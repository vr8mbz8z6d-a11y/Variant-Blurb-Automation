from __future__ import annotations
from models import VariantRecord
from variant_type import classify_variant_type
from myvariant_source import lookup_revel
from gnomad_source import lookup_gnomad
from clinvar_source import lookup_clinvar, find_co_located_variants
from spliceai_source import lookup_spliceai
from ensembl_hgvs_source import resolve_coordinates, resolve_coordinates_from_gene, get_splice_site_context
from nmd_predictor import predict_nmd
from hgvs_source_spare import HGVSResolutionError
from templates import build_blurb, normalize_clinical_classification


def annotate_variant(chrom: str, pos: int, ref: str, alt: str,
                      hgvs_c: str | None = None, hgvs_p: str | None = None,
                      gene: str | None = None, genome_build: str = "hg38",
                      hgvs_c_with_transcript: str | None = None,
                      clinical_classification: str | None = None,
                      debug: bool = False) -> VariantRecord:
    
    resolved_transcript: str | None = None
    resolved_strand: int | None = None

    if hgvs_c_with_transcript:
        try:
            resolved = resolve_coordinates(hgvs_c_with_transcript, genome_build=genome_build, debug=debug)
            chrom, pos, ref, alt = resolved.chrom, resolved.pos, resolved.ref, resolved.alt
            resolved_transcript = resolved.transcript
            resolved_strand = resolved.strand
            if debug:
                print(f"[pipeline] resolved {hgvs_c_with_transcript!r} -> "
                      f"{resolved.chrom}:{resolved.pos} {resolved.ref}>{resolved.alt}")
        except HGVSResolutionError as e:
            if debug:
                print(f"[pipeline] HGVS resolution failed, falling back to passed-in "
                      f"coordinates ({chrom}:{pos} {ref}>{alt}): {e}")

    record = VariantRecord(
        chrom=chrom, pos=pos, ref=ref, alt=alt,
        genome_build=genome_build, hgvs_c=hgvs_c, hgvs_p=hgvs_p, gene=gene,
        hgvs_transcript=resolved_transcript, gene_strand=resolved_strand,
        clinical_classification=clinical_classification,
    )
    record.variant_type = classify_variant_type(hgvs_p=hgvs_p, hgvs_c=hgvs_c)

    # gnomAD 
    record.gnomad = lookup_gnomad(chrom, pos, ref, alt, genome_build=genome_build, debug=debug)

    # REVEL 
    record.revel = lookup_revel(record.myvariant_hgvs_id(), genome_build=genome_build, debug=debug)

    # ClinVar&pmids
    ncbi_build = "GRCh38" if genome_build.lower() in ("hg38", "grch38") else "GRCh37"
    record.clinvar, record.pubmed_ids = lookup_clinvar(
    chrom,
    pos,
    ref,
    alt,
    genome_build=ncbi_build,
    hgvs=hgvs_c_with_transcript,
    debug=debug,
)

    # Co-located variants: OTHER ClinVar entries at the SAME amino acid
    # position as this one (e.g. query is p.Asn1355Lysfs*10, this finds
    # a separately-reported p.Asn1355Lys at the same residue) -- see
    # find_co_located_variants's docstring for the confirmed reference
    # example this reproduces. Needs both a transcript accession and a
    # protein-level hgvs_p to extract a search position from; silently
    # produces an empty list (no sentence rendered) if either is missing,
    # e.g. for a splice variant with no p. notation at all.
    transcript_for_colocated = resolved_transcript or (hgvs_c_with_transcript.split(":")[0]
                                                        if hgvs_c_with_transcript else None)
    if transcript_for_colocated and hgvs_p:
        record.co_located_variants = find_co_located_variants(
            transcript_for_colocated, hgvs_p,
            exclude_variation_id=record.clinvar.variation_id,
            debug=debug,
        )

    # SpliceAI
    # splice site variant.
    if record.variant_type in ("splice", "splice site variant"):
        hg_param = "38" if genome_build.lower() in ("hg38", "grch38") else "37"
        record.spliceai = lookup_spliceai(chrom, pos, ref, alt, genome_build=hg_param, debug=debug)

        # Splice-site exon/intron boundary context (e.g. "immediately
        # following exon 12"), for the descriptive sentence. Reuses the
        # transcript-qualified HGVS string used for coordinate
        # resolution; if that wasn't passed directly (e.g. the
        # gene+hgvs_c entry point), rebuild it from the resolved
        # transcript + hgvs_c.
        if record.variant_type == "splice":
            transcript_qualified = hgvs_c_with_transcript
            if not transcript_qualified and resolved_transcript and hgvs_c:
                transcript_qualified = f"{resolved_transcript}:{hgvs_c}"
            if transcript_qualified:
                record.splice_context = get_splice_site_context(
                    transcript_qualified, genome_build=genome_build, debug=debug
                )
            elif debug:
                print("[pipeline] no transcript-qualified HGVS available -- skipping "
                      "splice-site context lookup")

    # NMD (nonsense-mediated decay) prediction -- only for nonsense/
    # frameshift variants, same "build the transcript-qualified HGVS if
    # not already available" pattern as the splice-context lookup above.
    if record.variant_type in ("nonsense", "frameshift"):
        transcript_qualified = hgvs_c_with_transcript
        if not transcript_qualified and resolved_transcript and hgvs_c:
            transcript_qualified = f"{resolved_transcript}:{hgvs_c}"
        if transcript_qualified:
            record.nmd = predict_nmd(
                transcript_qualified, record.variant_type, hgvs_p=hgvs_p,
                genome_build=genome_build, debug=debug,
            )
        elif debug:
            print("[pipeline] no transcript-qualified HGVS available -- skipping "
                  "NMD prediction")

    return record


def annotate_and_render(chrom: str, pos: int, ref: str, alt: str,
                         hgvs_c: str | None = None, hgvs_p: str | None = None,
                         gene: str | None = None, genome_build: str = "hg38",
                         hgvs_c_with_transcript: str | None = None,
                         clinical_classification: str | None = None,
                         debug: bool = False) -> tuple[VariantRecord, str]:
    record = annotate_variant(chrom, pos, ref, alt, hgvs_c, hgvs_p, gene, genome_build,
                               hgvs_c_with_transcript=hgvs_c_with_transcript,
                               clinical_classification=clinical_classification,
                               debug=debug)
    return record, build_blurb(record)


def annotate_from_gene_and_hgvs_c(gene: str, hgvs_c: str, hgvs_p: str | None = None,
                                   genome_build: str = "hg38",
                                   clinical_classification: str | None = None,
                                   debug: bool = False) -> tuple[VariantRecord, str]:
    
    resolved = resolve_coordinates_from_gene(gene, hgvs_c, genome_build=genome_build, debug=debug)
    if debug:
        print(f"[pipeline] {gene} {hgvs_c} -> transcript {resolved.transcript} "
              f"({resolved.chrom}:{resolved.pos} {resolved.ref}>{resolved.alt})")

    record = annotate_variant(
        resolved.chrom, resolved.pos, resolved.ref, resolved.alt,
        hgvs_c=hgvs_c, hgvs_p=hgvs_p, gene=gene, genome_build=genome_build,
        clinical_classification=clinical_classification,
        debug=debug,
    )
    record.hgvs_transcript = resolved.transcript
    record.gene_strand = resolved.strand
    return record, build_blurb(record)


if __name__ == "__main__":
    from variant_input_parser import parse_combined_variant_string, VariantStringParseError

    print("Enter the variant below: ")
    print("(transcript accession is required; gene symbol and p. notation are optional)")
    raw_input_str = input("Variant: ").strip()

    try:
        parsed = parse_combined_variant_string(raw_input_str)
    except VariantStringParseError as e:
        print()
        print(f"Could not parse that input: {e}")
    else:
        genome_build = input("Genome build [hg38]: ").strip() or "hg38"
        clinical_classification_input = input(
            "Variant classification: "
            
        ).strip()
        clinical_classification = normalize_clinical_classification(clinical_classification_input)
        if clinical_classification_input and clinical_classification is None:
            print(f"Unrecognized variant classification {clinical_classification_input!r}; using ClinVar/automatic classification.")
        #debug_input = input("Show debug output for HGVS resolution? [y/N]: ").strip().lower()
        debug = False  #debug_input == "y"

        print()
        print(f"Parsed -- transcript: {parsed.transcript}, gene: {parsed.gene or '(not given)'}, "
              f"hgvs_c: {parsed.hgvs_c}, hgvs_p: {parsed.hgvs_p or '(not given)'}")

        try:
            
            if debug:
                print(f"[pipeline] resolving {parsed.hgvs_c_with_transcript!r} before "
                      f"running any downstream lookups...")
            resolved = resolve_coordinates(parsed.hgvs_c_with_transcript,
                                            genome_build=genome_build, debug=debug)
            record, blurb = annotate_and_render(
                chrom=resolved.chrom, pos=resolved.pos, ref=resolved.ref, alt=resolved.alt,
                hgvs_c=parsed.hgvs_c, hgvs_p=parsed.hgvs_p, gene=parsed.gene,
                genome_build=genome_build,
                hgvs_c_with_transcript=parsed.hgvs_c_with_transcript,
                clinical_classification=clinical_classification,
                debug=debug,
            )
        except HGVSResolutionError as e:
            print()
            print("Could not resolve genomic coordinates for this variant:")
            print(f"  {e}")
            print()
            print("Stopping here rather than guessing -- nothing downstream")
            print("(gnomAD/ClinVar/REVEL) was queried, so you won't get a confusing")
            print("blurb built from the wrong coordinates. Common causes:")
            print("  - malformed HGVS notation (e.g. missing the '>alt' half of a")
            print("    substitution, like 'c.1408+2T' instead of 'c.1408+2T>C')")
            print("  - a network timeout (rare; try again)")
            print("  - a transcript accession that doesn't match the gene")
            print("Re-run with debug=y to see the raw request/response if this persists.")
        else:
            print()
            print("Variant type:", record.variant_type)
            print("Resolved transcript:", record.hgvs_transcript)
            print()
            print("Blurb:")
            print(blurb if blurb else "(no facts resolved -- check network access / coordinates)")
