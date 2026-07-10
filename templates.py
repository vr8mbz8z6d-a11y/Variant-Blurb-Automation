from __future__ import annotations
import re
from models import VariantRecord


def gnomad_sentence(v: VariantRecord) -> str | None:
    g = v.gnomad
    if g.af is None:
        return None
    if (v.variant_type == "splice site variant" and g.af > 0
            and not v.clinvar.condition and not v.pubmed_ids):
        return None
    if g.af == 0:
        return "It was absent from large population studies such as the Genome Aggregation Database (gnomAD)."
    # Reference format: "present in 0.0001% (2/1,613,808; 0 homozygotes) of
    # total alleles in gnomAD". Thousands separators on the counts; hom
    # count always stated (0 is meaningful in a clinical read).
    count_str = ""
    if g.allele_count is not None and g.allele_number:
        count_str = f" ({g.allele_count:,}/{g.allele_number:,}; {g.homozygote_count or 0} homozygotes)"
    return (f"It is present in {g.af:.4%}{count_str} of total alleles in the "
            f"Genome Aggregation Database (gnomAD).")


def revel_sentence(v: VariantRecord) -> str | None:
    if v.variant_type == "splice site variant":
        return None
    if v.revel.score is None:
        return None
    return f"Computational prediction tools suggest an impact to protein function (REVEL {v.revel.score:.2f})."


def spliceai_sentence(v: VariantRecord) -> str | None:
    # Hard gate: only render for splice-classified variants, even if a
    # score is somehow present (e.g. from a stale/incorrect lookup).
    if v.variant_type != "splice":
        return None
    if v.spliceai.max_delta is None:
        return None
    return f"SpliceAI predicts a splicing impact with a maximum delta score of {v.spliceai.max_delta:.2f}."


_SPLICE_OFFSET_PATTERN = re.compile(r"(?P<sign>[+-])(?P<offset>\d+)[ACGT]>[ACGT]")


def splice_site_consensus_and_tools_sentence(v: VariantRecord) -> str | None:
    
    if v.variant_type != "splice site variant" or not v.hgvs_c:
        return None
    match = _SPLICE_OFFSET_PATTERN.search(v.hgvs_c)
    if not match:
        return None
    if v.gene_strand == 1:
        side = "5'"
    elif v.gene_strand == -1:
        side = "3'"
    else:
        return None
    consensus_clause = f"This variant is located in the {side} splice consensus region"
    if v.spliceai.max_delta is None:
        return f"{consensus_clause}."
    if v.spliceai.max_delta >= 0.2:
        tools_clause = "computational tools suggest an impact to splicing"
    else:
        tools_clause = "computational tools do not suggest an impact to splicing"
    return f"{consensus_clause}, and {tools_clause}."


def splice_context_sentence(v: VariantRecord) -> str | None:
    
    if v.variant_type != "splice":
        return None
    ctx = v.splice_context
    if ctx.offset is None or ctx.site_type is None or ctx.nearest_exon is None:
        return None
    relation = "following" if ctx.site_type == "donor" else "preceding"
    return (f"This variant alters the invariant {ctx.offset} nucleotide of the splice "
            f"{ctx.site_type} site immediately {relation} exon {ctx.nearest_exon}.")


def splice_mechanism_sentence(v: VariantRecord) -> str | None:
    
    if v.variant_type != "splice":
        return None
    ctx = v.splice_context
    if ctx.offset is None or ctx.site_type is None:
        return None
    if not v.gene:
        return None
    try:
        offset_magnitude = abs(int(ctx.offset))
    except (TypeError, ValueError):
        return None
    if offset_magnitude not in (1, 2):
        return None
    return (f"The \u00b11/\u00b12 {ctx.site_type} positions are essential for spliceosome "
            f"recognition, and disruption of this site is predicted to abolish normal "
            f"splicing, likely resulting in an abnormal or absent {v.gene} protein.")


_FRAMESHIFT_PATTERN = re.compile(
    r"p\.\(?[A-Za-z]{3}(?P<position>\d+)[A-Za-z]*fs(?:Ter|\*)(?P<stop_offset>\d+)"
)


def frameshift_sentence(v: VariantRecord) -> str | None:
    
    if v.variant_type != "frameshift":
        return None
    if not v.hgvs_p:
        return None
    match = _FRAMESHIFT_PATTERN.search(v.hgvs_p)
    if not match:
        return None
    position = match.group("position")
    stop_offset = match.group("stop_offset")
    return (f"This frameshift variant is predicted to alter the protein's amino acid "
            f"sequence beginning at position {position} and lead to a premature "
            f"termination codon {stop_offset} amino acids downstream. This alteration "
            f"is then predicted to lead to a truncated or absent protein.")


def nmd_sentence(v: VariantRecord) -> str | None:
    
    if v.variant_type not in ("nonsense", "frameshift"):
        return None
    if v.nmd.escapes_nmd is not True:
        # Covers both escapes_nmd is None (unresolved -- don't guess) and
        # escapes_nmd is False (triggers NMD -- say nothing, per request).
        return None
    return (f"This variant is predicted to escape nonsense-mediated decay (NMD), "
            f"as {v.nmd.reason}, and is therefore more likely to result in a "
            f"truncated protein product rather than complete loss of protein "
            f"expression.")


def clinvar_other_labs_sentence(v: VariantRecord) -> str | None:
    

    # CASE 1: Truly no ClinVar record
    if not v.clinvar.variation_id:
        return "This variant has not been reported in ClinVar."

    id_str = f" (ClinVar ID: {v.clinvar.variation_id})"
    fallback_sentence = (
        f"This variant is present in ClinVar{id_str}, but has no "
        f"submitted clinical interpretations."
    )

    # CASE 2: Record exists but no meaningful data at all
    if not v.clinvar.submissions and not v.clinvar.aggregate_classification:
        return fallback_sentence

    # CASE 3: Build a label from the aggregate call, or from distinct
    
    label = (v.clinvar.aggregate_classification or "").strip()
    if not label:
        seen = []
        for s in v.clinvar.submissions:
            c = (s.classification or "").strip().capitalize()
            if c and c not in seen:  # skip blank/whitespace-only entries
                seen.append(c)
        label = "/ ".join(seen)

    # Safety net: if nothing usable was found (e.g. every submission's
    # classification text was blank)
    if not label:
        return fallback_sentence

    return f"This variant has been classified as {label} by other clinical laboratories{id_str}."


def co_located_variants_sentence(v: VariantRecord) -> str | None:
    """
    Reports OTHER ClinVar variants affecting the SAME amino acid
    position as this one, e.g.:
      "Another variant in the same position, p.Arg5105Gln, has been
       classified as Conflicting classifications of pathogenicity by
       other clinical laboratories (ClinVar ID: 69441)."

    Confirmed against a real reference report (wording: "Another variant
    in the same position, p.Arg5105Gln, has conflicting classifications
    of pathogenicity by other clinical laboratories (ClinVar ID:
    69441).") -- one deliberate wording choice worth flagging: this
    reuses the exact "has been classified as X by other clinical
    laboratories" phrase already established in
    clinvar_other_labs_sentence, rather than the reference's slightly
    more informal "has X" shortcut, since that shortcut only reads
    naturally for this one classification label ("conflicting
    classifications...") and wouldn't generalize cleanly to every
    possible label (e.g. "has Pathogenic" reads oddly; "has been
    classified as Pathogenic" doesn't).

    Renders ONE sentence per co-located variant found (up to the cap
    already applied in find_co_located_variants), so multiple co-located
    variants produce multiple consecutive sentences rather than one
    long, comma-stacked sentence -- no reference was given for the
    multiple-variant case, so this is a reasonable default, not a
    confirmed convention.

    Omitted entirely if none were found (the common case) -- this data
    is populated by find_co_located_variants in clinvar_source.py.
    """
    if not v.co_located_variants:
        return None

    sentences = []
    for cv in v.co_located_variants:
        id_str = f" (ClinVar ID: {cv.clinvar_id})" if cv.clinvar_id else ""
        if cv.classification:
            sentences.append(
                f"Another variant in the same position, {cv.protein_change}, has been "
                f"classified as {cv.classification} by other clinical laboratories{id_str}."
            )
        else:
            sentences.append(
                f"Another variant in the same position, {cv.protein_change}, is also "
                f"reported in ClinVar{id_str}."
            )
    return " ".join(sentences)


MAX_PMIDS_IN_BLURB = 3  # cap PMIDs shown in the blurb (you asked for 2-4)


def _protein_change_phrase(v: VariantRecord) -> str:
    
    if v.variant_type == "splice site variant":
        type_word = " splice site"
    else:
        type_word = f" {v.variant_type}" if v.variant_type and v.variant_type != "unknown" else ""
    gene_phrase = f" in {v.gene}" if v.gene else ""

    if v.hgvs_p:
        return f"The {v.hgvs_p}{type_word} variant{gene_phrase}"
    if v.hgvs_c:
        if v.variant_type == "splice site variant":
            return f"The {v.hgvs_c} splice site variant{gene_phrase}"
        # Coding notation alone doesn't carry a "missense/nonsense/splice"
        # type word naturally the way protein notation does (e.g. "The
        # c.1838+2T>C splice variant" reads a bit awkwardly stacked with
        # a c. string) -- the confirmed reference wording omits the type
        # word here, so we do too, e.g. "The c.1838+2T>C variant in P3H1".
        return f"The {v.hgvs_c} variant{gene_phrase}"
    # no HGVS notation available at all -> fall back to a generic subject
    return f"This{type_word} variant{gene_phrase}"


def opening_sentence(v: VariantRecord) -> str | None:
   
    subject = _protein_change_phrase(v)

    if (v.variant_type == "splice site variant" and not v.clinvar.condition
            and not v.pubmed_ids and v.gnomad.af is not None and v.gnomad.af > 0):
        count_str = ""
        if v.gnomad.allele_count is not None and v.gnomad.allele_number:
            count_str = (
                f" ({v.gnomad.allele_count:,}/{v.gnomad.allele_number:,}; "
                f"{v.gnomad.homozygote_count or 0} homozygotes)"
            )
        return (f"{subject} has not been previously reported in affected individuals "
                f"but was identified in {v.gnomad.af:.4%}{count_str} of total "
                f"alleles in the Genome Aggregation Database (gnomAD).")

    # condition clause
    if v.clinvar.condition:
        context = f"has previously been reported in individuals with {v.clinvar.condition}"
    elif v.pubmed_ids:
        context = "has previously been reported in the literature"
    else:
        return f"{subject} has not been previously reported in affected individuals."

    # PMID clause (capped, newest-first) -- reuses the same selection logic
    pmid_clause = ""
    if v.pubmed_ids:
        try:
            ordered = sorted(v.pubmed_ids, key=lambda x: int(x), reverse=True)
        except ValueError:
            ordered = list(v.pubmed_ids)
        shown = ordered[:MAX_PMIDS_IN_BLURB]
        suffix = " among others" if len(v.pubmed_ids) > len(shown) else ""
        pmid_clause = f" (PMID: {', '.join(shown)}{suffix})"

    return f"{subject} {context}{pmid_clause}."


def functional_study_sentence(v: VariantRecord) -> str | None:
    # Deliberately a placeholder, not a real claim -- see FunctionalStudyFlag
    # in models.py. This pipeline does not yet assert functional evidence.
    return None


# Substrings (checked case-insensitively) that indicate genuine clinical
# uncertainty rather than a confident benign/pathogenic call. Using
# substring matching, not exact equality, since ClinVar's own wording has
# varied over time -- e.g. it used to say "conflicting interpretations of
# pathogenicity" before switching to "conflicting classifications of
# pathogenicity". Both are caught by checking for "conflicting".
_UNCERTAIN_CLASSIFICATION_MARKERS = ("uncertain significance", "conflicting")


_CLINICAL_CLASSIFICATION_ALIASES = {
    "pathogenic": "Pathogenic",
    "pathogeneic": "Pathogenic",
    "likely pathogenic": "Likely pathogenic",
    "likely pathogeneic": "Likely pathogenic",
    "likely-pathogenic": "Likely pathogenic",
    "lp": "Likely pathogenic",
    "vus": "VUS",
    "uncertain significance": "VUS",
    "variant of uncertain significance": "VUS",
    "unknown significance": "VUS",
    "likely benign": "Likely benign",
    "likely-benign": "Likely benign",
    "lb": "Likely benign",
    "benign": "Benign",
    "conflicting": "Conflicting classifications of pathogenicity",
    "conflicting classifications": "Conflicting classifications of pathogenicity",
    "conflicting classifications of pathogenicity": "Conflicting classifications of pathogenicity",
}


def normalize_clinical_classification(value: str | None) -> str | None:
    if value is None:
        return None
    key = re.sub(r"\s+", " ", value.strip().lower())
    if not key:
        return None
    return _CLINICAL_CLASSIFICATION_ALIASES.get(key)


def closing_summary_sentence(v: VariantRecord) -> str | None:
    
    label = normalize_clinical_classification(v.clinical_classification)
    if label is None:
        label = v.clinvar.aggregate_classification
    if not label:
        return None

    label_lower = label.lower()
    if label == "VUS" or any(marker in label_lower for marker in _UNCERTAIN_CLASSIFICATION_MARKERS):
        return ("In summary, additional information is needed to fully assess "
                "the clinical significance of this variant.")

    return f"In summary, this variant meets our criteria to be classified as {label}."


# Order mirrors the reference report:
#   opening (identity + literature) -> gnomAD -> in-silico -> functional
#   (currently a no-op) -> ClinVar other labs -> closing summary
SENTENCE_MODULES = [
    opening_sentence,
    gnomad_sentence,
    revel_sentence,
    spliceai_sentence,
    splice_site_consensus_and_tools_sentence,
    splice_context_sentence,
    splice_mechanism_sentence,
    frameshift_sentence,
    nmd_sentence,
    functional_study_sentence,
    clinvar_other_labs_sentence,
    co_located_variants_sentence,
    closing_summary_sentence,
]


def build_blurb(v: VariantRecord) -> str:
    sentences = [fn(v) for fn in SENTENCE_MODULES]
    return " ".join(s for s in sentences if s)
