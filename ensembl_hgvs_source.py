"""
HGVS coding-variant -> genomic coordinate resolution via Ensembl's VEP
REST API, replacing VariantValidator (hgvs_source_spare.py) entirely
after repeated live timeouts on that service.

WHY ENSEMBL: rest.ensembl.org is core EMBL-EBI infrastructure, used by a
huge number of downstream tools; VariantValidator is a much smaller
academic tool from a single university group, and had both a documented
stale-cache bug and repeated live timeouts during this project's testing.
This is not a guarantee of higher uptime, but it's a reasonable basis for
trying it as the primary resolver.

VERIFIED (from Ensembl's own documentation and a real Bioconductor
package usage example -- not guessed):
  - Endpoint: GET https://rest.ensembl.org/vep/human/hgvs/{hgvs}
    (GRCh38). GRCh37 uses a separate host: https://grch37.rest.ensembl.org
    -- NOT a query parameter, an entirely different subdomain.
  - Explicitly supports RefSeq transcript accessions (NM_...) as input,
    not just Ensembl's own ENST... IDs -- confirmed via Ensembl's own
    documented examples (e.g. "NM_005239.6:c.190G>A").
  - Response is a JSON LIST (even for one input variant). Confirmed
    top-level fields on each element (from a real Bioconductor
    ensemblVEP vignette showing actual field names from a live call):
    seq_region_name, start, end, allele_string, strand, assembly_name,
    input, id, most_severe_consequence, transcript_consequences,
    colocated_variants, regulatory_feature_consequences.
  - allele_string is "REF/ALT" but is NOT reliably on the genomic
    (plus-strand, VCF-style) orientation. CORRECTED after a live bug: for
    a minus-strand gene, allele_string comes back in TRANSCRIPT-sense
    (cDNA) orientation, matching the raw c. notation, not the reference
    genome strand. Confirmed against a real case: ASPM
    NM_018136.5:c.8820+2T>C (a minus-strand gene) returned allele_string
    "T/C" with strand=-1, while the true plus-strand change -- confirmed
    independently via spliceailookup.broadinstitute.org's own resolution
    -- is A>G (T>C's exact reverse-complement). This module now reads the
    top-level `strand` field and reverse-complements ref/alt whenever
    strand == -1, before returning coordinates to any caller.

NOT YET LIVE-VERIFIED (flagging honestly rather than guessing further):
  - Exact behavior for intronic splice-site offset notation
    (e.g. c.1408+2T>C) specifically -- VEP is built to handle this, but
    no live example with this exact notation pattern was found during
    research. If this fails on splice variants, that's the first thing
    to check with debug=True.
  - Insertion coordinate convention: Ensembl's start/end for pure
    insertions can follow a different convention than VCF (their docs
    describe insertions as start > end, i.e. start = end + 1). This
    module currently just uses `start` as-is, which is correct for SNVs
    and substitutions (start == end) but has NOT been verified correct
    for pure insertions. If you hit an insertion variant and the
    resulting position looks off by one, this is the first place to look.
"""
from __future__ import annotations
from typing import Optional
import re
import time
import requests
from urllib.parse import quote

from hgvs_source_spare import HGVSResolutionError, ResolvedCoordinates

GRCH38_BASE_URL = "https://rest.ensembl.org"
GRCH37_BASE_URL = "https://grch37.rest.ensembl.org"
DEFAULT_TIMEOUT_SECONDS = 40

# Ensembl's public REST API has documented, recurring periods of slowness
# under server load (not specific to any one gene) -- see e.g.
# lists.ensembl.org/pipermail/dev_ensembl.org/2018-July/008022.html and
# github.com/Ensembl/ensembl-vep/issues/1807. Separately, and NOT
# confirmed against a written Ensembl source (flagging honestly): large,
# exon-dense genes like TTN (~364 exons, many overlapping transcripts)
# are plausibly slower for VEP to compute consequences across than a
# typical gene, since /vep/human/hgvs/ annotates every overlapping
# transcript. Both point the same direction -- a timeout here is more
# likely "needs more time" than "something is broken" -- so retry with a
# longer timeout each attempt rather than failing on the first one.
MAX_RETRIES = 2  # total attempts = MAX_RETRIES + 1
RETRY_BACKOFF_SECONDS = 3
RETRY_TIMEOUT_MULTIPLIER = 1.7  # give a slow gene more room on each retry

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _reverse_complement(seq: str) -> Optional[str]:
    """
    Reverse-complement a base string (A<->T, C<->G). Returns None if the
    sequence contains anything outside standard ACGT, since we should
    never silently guess a complement for something we don't recognize
    (e.g. indel placeholder characters) -- callers must treat None as
    "cannot safely strand-correct this one, don't try."
    """
    seq = seq.upper()
    if not seq or any(base not in _COMPLEMENT for base in seq):
        return None
    return "".join(_COMPLEMENT[base] for base in reversed(seq))


def _parse_transcript_and_hgvs(hgvs_c_with_transcript: str) -> tuple[str, str]:
    """Split 'NM_000091.5:c.1408+2T>C' into ('NM_000091.5', 'c.1408+2T>C').
    Requires an explicit transcript accession -- never guesses one."""
    if ":" not in hgvs_c_with_transcript:
        raise HGVSResolutionError(
            f"Expected a transcript-qualified HGVS string like "
            f"'NM_000162.5:c.1021G>A', got {hgvs_c_with_transcript!r} "
            f"(no transcript accession found before the colon)."
        )
    transcript, hgvs_c = hgvs_c_with_transcript.split(":", 1)
    transcript = transcript.strip()
    hgvs_c = hgvs_c.strip()
    if not transcript or not hgvs_c.startswith("c."):
        raise HGVSResolutionError(
            f"Could not parse {hgvs_c_with_transcript!r} as "
            f"'TRANSCRIPT:c.change'."
        )
    return transcript, hgvs_c


def _base_url_for_build(genome_build: str) -> str:
    if genome_build.lower() in ("hg19", "grch37", "37"):
        return GRCH37_BASE_URL
    return GRCH38_BASE_URL


def _verify_response_matches_query(result: dict, hgvs_c_with_transcript: str,
                                    debug: bool = False) -> None:
    """
    Stale/wrong-match defense, mirroring the same discipline used
    elsewhere in this pipeline (hgvs_source_spare.py's _verify_match, and
    the ClinVar HGVS-search fix) after finding a real case where an
    unverified text-based match returned data for the WRONG variant.
    Ensembl's response echoes the exact input string back in the "input"
    field -- if it doesn't match what we sent, something is wrong and we
    must not trust the coordinates.
    """
    echoed = result.get("input")
    if echoed is None:
        raise HGVSResolutionError(
            f"Ensembl response has no 'input' field to verify against -- "
            f"refusing to trust unverified coordinates. Full result keys: "
            f"{list(result.keys())!r}"
        )
    # Normalize whitespace/case for comparison; Ensembl echoes the input
    # essentially verbatim, so this should be a near-exact match.
    if echoed.strip().lower() != hgvs_c_with_transcript.strip().lower():
        raise HGVSResolutionError(
            f"Ensembl response's echoed input ({echoed!r}) does not match "
            f"the query ({hgvs_c_with_transcript!r}) -- refusing to trust "
            f"this result rather than risk returning data for the wrong "
            f"variant."
        )
    if debug:
        print(f"[ensembl_hgvs_source] verified: response echoes back the exact query")


def _get_with_retries(http: "requests.Session", url: str, params: dict,
                       timeout_seconds: float, debug: bool = False) -> "requests.Response":
    """
    GET with retries ONLY on timeout/connection errors -- NOT on a real
    HTTP error response (4xx/5xx), where retrying is pointless (a
    rejected/malformed query won't become valid on attempt 2). Timeouts
    and connection errors, by contrast, are the documented signature of
    Ensembl's public REST API under transient load (see the module
    docstring), so retrying with a longer timeout is worth doing before
    giving up.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, MAX_RETRIES + 2):
        try:
            if debug:
                print(f"[ensembl_hgvs_source] GET attempt {attempt}/{MAX_RETRIES + 1} "
                      f"(timeout={timeout_seconds:.0f}s): {url}")
            return http.get(url, params=params, timeout=timeout_seconds)
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_exc = exc
            if debug:
                print(f"[ensembl_hgvs_source] attempt {attempt} failed: {exc}")
            if attempt <= MAX_RETRIES:
                if debug:
                    print(f"[ensembl_hgvs_source] retrying in {RETRY_BACKOFF_SECONDS}s "
                          f"with a longer timeout...")
                time.sleep(RETRY_BACKOFF_SECONDS)
                timeout_seconds *= RETRY_TIMEOUT_MULTIPLIER
    raise last_exc


def _fetch_reference_base(chrom: str, pos: int, genome_build: str = "GRCh38",
                           timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
                           session: Optional["requests.Session"] = None,
                           debug: bool = False) -> Optional[str]:
    """
    Fetch a single reference base at a genomic position via Ensembl's
    sequence endpoint. Used to construct a proper VCF-style anchor base
    for indel normalization (see _normalize_indel_to_vcf).

    CONFIRMED endpoint and region format from Ensembl's own documentation
    (rest.ensembl.org/documentation/info/sequence_region):
        GET /sequence/region/human/{chrom}:{start}..{end}

    NOT YET LIVE-VERIFIED (flagging honestly): the exact JSON field name
    for the returned sequence when requesting content-type=application/
    json. Assumed to be {"seq": "..."}, the standard Ensembl REST
    convention for sequence-returning endpoints, but not confirmed
    against a live response for this specific endpoint. If this returns
    None unexpectedly, run with debug=True and check the raw response
    printed -- the field name is the first thing to adjust.
    """
    base_url = _base_url_for_build(genome_build)
    chrom_clean = chrom.replace("chr", "").replace("Chr", "")
    region = f"{chrom_clean}:{pos}..{pos}"
    url = f"{base_url}/sequence/region/human/{region}"
    params = {"content-type": "application/json"}
    http = session or requests.Session()

    if debug:
        print(f"[ensembl_hgvs_source] GET {url} params={params} "
              f"(fetching anchor base for indel normalization)")

    try:
        response = _get_with_retries(http, url, params, timeout_seconds, debug=debug)
    except requests.RequestException as exc:
        if debug:
            print(f"[ensembl_hgvs_source] reference base fetch failed: {exc}")
        return None

    if debug:
        print(f"[ensembl_hgvs_source] reference base fetch status {response.status_code}")
        print(f"[ensembl_hgvs_source] reference base fetch response: {response.text[:300]}")

    if response.status_code != 200:
        return None

    try:
        payload = response.json()
        seq = payload.get("seq") if isinstance(payload, dict) else None
    except ValueError:
        return None

    if not seq or len(seq) != 1:
        if debug:
            print(f"[ensembl_hgvs_source] unexpected sequence response shape "
                  f"(expected a single-base 'seq' field): {payload!r}")
        return None

    return seq.upper()


def _normalize_indel_to_vcf(chrom: str, pos: int, ref: str, alt: str,
                             genome_build: str = "GRCh38",
                             timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
                             session: Optional["requests.Session"] = None,
                             debug: bool = False) -> tuple[str, int, str, str]:
    """
    Convert Ensembl's indel convention (a literal "-" for the absent
    allele, e.g. allele_string "A/-" for a deletion) into standard
    VCF-style anchor-base representation, which is what every downstream
    consumer (gnomAD, MyVariant.info, ClinVar's coordinate search) needs.

    CONFIRMED real bug this fixes: querying NF1 c.5920del produced
    allele_string "A/-". Passing ref="A", alt="-" straight through built
    a malformed gnomAD query "17-31334945-A--" (a trailing double-
    hyphen), which gnomAD correctly rejected with "Invalid variant ID"
    -- silently and permanently blocking the "absent from gnomAD"
    sentence for every deletion, and separately causing MyVariant.info
    to 404 on the same malformed representation.

    Standard VCF convention:
      Deletion of base(s) at position P: anchor = reference base at P-1
        (unaffected by the deletion). VCF: POS=P-1, REF=anchor+deleted,
        ALT=anchor.
      Insertion (Ensembl represents this with start = end+1): anchor =
        reference base at (start-1). VCF: POS=start-1, REF=anchor,
        ALT=anchor+inserted.

    SCOPE: works for BOTH plus- and minus-strand indels. This is called
    (see resolve_coordinates) AFTER strand correction, which now
    correctly reverse-complements only the real-sequence side of a
    minus-strand indel (leaving the "-" placeholder for "absence"
    unchanged, since absence has no meaningful complement). Once strand
    correction has run, ref/alt arriving here are always in genomic
    plus-strand sense, and the anchor-base position formula below
    (pos-1) is correct regardless of the gene's strand -- it's a
    genomic-adjacency concept, not a transcript-relative one, since
    Ensembl's start/end are always given in genomic coordinates. (An
    earlier version of this comment incorrectly claimed minus-strand
    indels would need the anchor fetched from pos+1 instead -- that was
    an unverified guess made before working through the actual geometry
    with a real case; it was wrong, and has been corrected here.)

    Returns (chrom, pos, ref, alt) unchanged if this isn't an indel
    (neither ref nor alt is "-").

    Raises:
        HGVSResolutionError: if the anchor base can't be fetched --
            refuses to guess or pass through a malformed ID rather than
            send downstream lookups a broken query.
    """
    if ref != "-" and alt != "-":
        return chrom, pos, ref, alt  # not an indel, nothing to normalize

    anchor_pos = pos - 1
    anchor_base = _fetch_reference_base(chrom, anchor_pos, genome_build=genome_build,
                                         timeout_seconds=timeout_seconds, session=session, debug=debug)
    if anchor_base is None:
        raise HGVSResolutionError(
            f"Could not fetch the reference anchor base at {chrom}:{anchor_pos} "
            f"needed to build a proper VCF-style representation of this indel "
            f"(Ensembl gave ref={ref!r}, alt={alt!r} at {chrom}:{pos}) -- refusing "
            f"to guess a malformed variant ID rather than send gnomAD/ClinVar/"
            f"MyVariant.info a broken query."
        )

    if alt == "-":
        # Deletion: ref is the deleted sequence itself (from allele_string).
        new_ref = anchor_base + ref
        new_alt = anchor_base
    else:
        # Insertion: alt is the inserted sequence itself.
        new_ref = anchor_base
        new_alt = anchor_base + alt

    if debug:
        print(f"[ensembl_hgvs_source] normalized indel to VCF style: "
              f"{chrom}:{pos} {ref}>{alt} -> {chrom}:{anchor_pos} {new_ref}>{new_alt}")

    return chrom, anchor_pos, new_ref, new_alt


def resolve_coordinates(
    hgvs_c_with_transcript: str,
    genome_build: str = "GRCh38",
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> ResolvedCoordinates:
    """
    Resolve a transcript-qualified HGVS coding variant description to
    strand-corrected genomic coordinates via Ensembl's VEP REST API.

    Drop-in replacement for hgvs_source_spare.resolve_coordinates() --
    same signature shape, same return type, same exception type -- so
    callers (pipeline.py) don't need structural changes, just a
    different import.

    Raises:
        HGVSResolutionError: on malformed input, network/HTTP failure,
            unexpected response shape, or a response that fails the
            input-echo verification check.
    """
    transcript, hgvs_c = _parse_transcript_and_hgvs(hgvs_c_with_transcript)
    base_url = _base_url_for_build(genome_build)
    http = session or requests.Session()

    url = f"{base_url}/vep/human/hgvs/{quote(hgvs_c_with_transcript, safe='')}"
    params = {"content-type": "application/json"}

    if debug:
        print(f"[ensembl_hgvs_source] GET {url} params={params}")

    try:
        response = _get_with_retries(http, url, params, timeout_seconds, debug=debug)
    except requests.Timeout as exc:
        raise HGVSResolutionError(
            f"Ensembl VEP timed out for {hgvs_c_with_transcript!r} after "
            f"{MAX_RETRIES + 1} attempt(s) (final timeout "
            f"{timeout_seconds * (RETRY_TIMEOUT_MULTIPLIER ** MAX_RETRIES):.0f}s): {exc}. "
            f"Ensembl's public REST API has documented periods of slowness under "
            f"load, and large, exon-dense genes (TTN especially) plausibly take "
            f"VEP longer to annotate than a typical gene. If this keeps happening, "
            f"try again in a few minutes, or pass a longer timeout_seconds to "
            f"resolve_coordinates() directly."
        ) from exc
    except requests.RequestException as exc:
        raise HGVSResolutionError(
            f"Network error calling Ensembl VEP for {hgvs_c_with_transcript!r}: {exc}"
        ) from exc

    if debug:
        print(f"[ensembl_hgvs_source] status {response.status_code}")
        print(f"[ensembl_hgvs_source] response (first 600 chars): {response.text[:600]}")

    if response.status_code != 200:
        raise HGVSResolutionError(
            f"Ensembl VEP returned HTTP {response.status_code} for "
            f"{hgvs_c_with_transcript!r}: {response.text[:500]}"
        )

    try:
        payload = response.json()
    except ValueError as exc:
        raise HGVSResolutionError(
            f"Ensembl VEP response was not valid JSON for {hgvs_c_with_transcript!r}"
        ) from exc

    # Ensembl returns a dict with an "error" key on failure, and a LIST of
    # result objects on success (even for a single input variant).
    if isinstance(payload, dict) and "error" in payload:
        raise HGVSResolutionError(
            f"Ensembl VEP rejected {hgvs_c_with_transcript!r}: {payload['error']}"
        )
    if not isinstance(payload, list) or not payload:
        raise HGVSResolutionError(
            f"Ensembl VEP returned an unexpected response shape for "
            f"{hgvs_c_with_transcript!r}: {payload!r}"
        )
    if len(payload) > 1 and debug:
        print(f"[ensembl_hgvs_source] WARNING: expected 1 result, got "
              f"{len(payload)} -- using the first one")

    result = payload[0]
    _verify_response_matches_query(result, hgvs_c_with_transcript, debug=debug)

    try:
        chrom = str(result["seq_region_name"])
        pos = int(result["start"])
        allele_string = result["allele_string"]
        strand = result.get("strand")
    except (KeyError, ValueError, TypeError) as exc:
        raise HGVSResolutionError(
            f"Ensembl VEP result is missing expected fields for "
            f"{hgvs_c_with_transcript!r}: {result!r}"
        ) from exc

    if "/" not in allele_string:
        raise HGVSResolutionError(
            f"Ensembl VEP 'allele_string' was not in the expected "
            f"'REF/ALT' format: {allele_string!r}"
        )
    ref, alt = allele_string.split("/", 1)

    # CONFIRMED BUG FIX: Ensembl's allele_string is NOT always already on
    # the plus-strand/VCF orientation as originally assumed -- for a
    # minus-strand gene (strand == -1), it stays in transcript-sense
    # (cDNA) orientation. Verified against a live case: ASPM
    # NM_018136.5:c.8820+2T>C returned allele_string "T/C" with
    # strand=-1, but the true plus-strand genomic change (confirmed
    # independently via spliceailookup.broadinstitute.org's own
    # resolution) is A>G -- the exact reverse-complement of T>C. Every
    # downstream lookup (gnomAD, MyVariant, SpliceAI) needs the
    # plus-strand orientation, so we correct it here rather than pass
    # transcript-sense alleles downstream.
    #
    # INDEL FIX: for a deletion/insertion, one side of allele_string is
    # Ensembl's literal "-" placeholder for "nothing" (e.g. "TCAA/-" for
    # a deletion). "-" represents absence, which has no meaningful
    # complement and is the SAME regardless of strand -- it should pass
    # through unchanged, not be rejected. Only the real sequence side
    # needs reverse-complementing. Confirmed correct by working through
    # the actual geometry with real numbers (BRCA1, a minus-strand gene,
    # c.4065_4068delTCAA): the genomic-order deleted sequence is exactly
    # reverse_complement("TCAA") = "TTGA", and this plugs directly into
    # the existing plus-strand indel-normalization logic (position math
    # is unaffected -- an EARLIER version of this comment incorrectly
    # assumed the anchor-base position formula also needed to change for
    # minus-strand indels; it doesn't, because VCF anchor-base placement
    # is a genomic-adjacency concept, and Ensembl's start/end are already
    # given in genomic (not transcript-relative) terms regardless of
    # gene strand).
    if strand in (-1, "-1"):
        rc_ref = ref if ref == "-" else _reverse_complement(ref)
        rc_alt = alt if alt == "-" else _reverse_complement(alt)
        if rc_ref is None or rc_alt is None:
            raise HGVSResolutionError(
                f"Variant is on the minus strand (strand=-1) but ref/alt "
                f"{ref!r}/{alt!r} contain characters this module can't "
                f"safely reverse-complement. Refusing to guess -- verify "
                f"this variant's genomic orientation manually before "
                f"trusting any downstream lookup for it."
            )
        if debug:
            print(f"[ensembl_hgvs_source] minus-strand gene (strand=-1): "
                  f"correcting allele_string {ref}>{alt} (transcript-sense) "
                  f"to {rc_ref}>{rc_alt} (plus-strand/genomic)")
        ref, alt = rc_ref, rc_alt
    elif debug:
        print(f"[ensembl_hgvs_source] plus-strand gene (strand={strand!r}): "
              f"allele_string {ref}>{alt} used as-is")

    # Indel normalization: Ensembl represents a deletion/insertion's
    # absent allele with a literal "-" (e.g. "A/-"), which is NOT valid
    # VCF-style representation -- gnomAD, MyVariant.info, and ClinVar's
    # coordinate search all expect a proper anchor-base representation
    # instead. See _normalize_indel_to_vcf's docstring for the confirmed
    # bug this fixes. No-op for plain substitutions (neither ref nor alt
    # is "-").
    chrom, pos, ref, alt = _normalize_indel_to_vcf(
        chrom, pos, ref, alt, genome_build=genome_build,
        timeout_seconds=timeout_seconds, session=http, debug=debug,
    )

    if debug:
        print(f"[ensembl_hgvs_source] resolved: {chrom}:{pos} {ref}>{alt} "
              f"(from allele_string {allele_string!r}, strand={strand!r})")

    return ResolvedCoordinates(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        genome_build=genome_build,
        transcript=transcript,
        hgvs_c=hgvs_c,
        validated_description=result.get("input", hgvs_c_with_transcript),
        strand=int(strand) if strand in (-1, 1, "-1", "1") else None,
    )


def resolve_coordinates_from_gene(
    gene_symbol: str,
    hgvs_c: str,
    genome_build: str = "GRCh38",
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> ResolvedCoordinates:
    """
    LOWER CONFIDENCE THAN resolve_coordinates(). Resolves a gene symbol's
    Ensembl Canonical transcript via the /lookup/symbol endpoint and uses
    that as the transcript for coordinate resolution.

    IMPORTANT CAVEAT: Ensembl Canonical is NOT guaranteed identical to
    the RefSeq MANE Select transcript ClinVar/clinical labs typically
    use. Per the MANE project's own published analysis (Morales et al.
    2022, Nature), the transcript ClinVar submitters actually used
    differed from MANE Select for about 7% of reviewed genes, and
    gnomAD's canonical differed for about 19%. For clinical use, prefer
    resolve_coordinates() with an explicit, known transcript accession
    (e.g. from a VarSeq export) whenever you have one -- this function is
    a convenience fallback, not the recommended default.
    """
    http = session or requests.Session()
    base_url = _base_url_for_build(genome_build)
    url = f"{base_url}/lookup/symbol/homo_sapiens/{quote(gene_symbol, safe='')}"
    params = {"content-type": "application/json", "expand": 0}

    if debug:
        print(f"[ensembl_hgvs_source] GET {url} params={params}")

    try:
        response = http.get(url, params=params, timeout=timeout_seconds)
    except requests.RequestException as exc:
        raise HGVSResolutionError(
            f"Network error calling Ensembl lookup/symbol for {gene_symbol!r}: {exc}"
        ) from exc

    if response.status_code != 200:
        raise HGVSResolutionError(
            f"Ensembl lookup/symbol returned HTTP {response.status_code} for "
            f"{gene_symbol!r}: {response.text[:500]}"
        )

    try:
        gene_data = response.json()
    except ValueError as exc:
        raise HGVSResolutionError(
            f"Ensembl lookup/symbol response was not valid JSON for {gene_symbol!r}"
        ) from exc

    transcript = gene_data.get("canonical_transcript")
    if not transcript:
        raise HGVSResolutionError(
            f"No canonical_transcript found for gene {gene_symbol!r} in "
            f"Ensembl's response: {gene_data!r}"
        )
    # canonical_transcript often includes a trailing Ensembl-internal
    # version suffix like "ENST00000544455.6.1" -- strip it back to the
    # standard "ENST00000544455.6" form.
    transcript = re.sub(r"^(ENST\d+\.\d+)\.\d+$", r"\1", transcript)

    if debug:
        print(f"[ensembl_hgvs_source] Ensembl Canonical transcript for "
              f"{gene_symbol!r}: {transcript} (NOT guaranteed == MANE Select)")

    return resolve_coordinates(
        f"{transcript}:{hgvs_c.strip()}",
        genome_build=genome_build,
        timeout_seconds=timeout_seconds,
        session=session,
        debug=debug,
    )


# ---------------------------------------------------------------------------
# Splice-site exon/intron boundary lookup
#
# Purpose: build a sentence like "This variant alters the invariant +2
# nucleotide of the splice donor site immediately following exon 12."
# for splice-region variants (c.NNNN+k or c.NNNN-k notation).
#
# CONFIRMED (from Ensembl's own official VEP documentation, not a guess):
# VEP reports which exon/intron a variant falls in as a "N/Total" string
# (e.g. "12/38" = intron 12 of 38), available via the `numbers` option
# ("Equivalent to --numbers" in Ensembl's docs).
#
# NOT YET LIVE-VERIFIED (flagging honestly): the exact REST query
# parameter names/values needed to turn this on for the
# /vep/human/hgvs/... endpoint specifically. This module passes
# `numbers=1` (to get the exon/intron field at all) and `hgvs=1` (to get
# a per-transcript `hgvsc` string back, needed to pick the correct
# transcript_consequences entry when VEP reports consequences across
# several overlapping transcripts of the same gene -- confirmed this
# happens in practice, e.g. ASPM returned consequences for BOTH
# ENST00000294732 (protein-coding, the one actually queried) and
# ENST00000367408 (a non-coding transcript) for the same HGVS query).
# If this comes back without exon/intron data, run with debug=True and
# check the raw payload -- the parameter names are the first thing to
# adjust.
# ---------------------------------------------------------------------------

def _parse_hgvs_intronic_offset(hgvs_c: str) -> Optional[str]:
    """
    Extract the intronic offset from HGVS c. notation, e.g.
    'c.1838+2T>C' -> '+2', 'c.1839-2A>G' -> '-2'. Returns None for
    non-intronic (exonic) HGVS descriptions -- this function is only
    meaningful for splice-region variants.
    """
    match = re.search(r"[+-]\d+(?=[ACGTacgt]>|del|ins|dup)", hgvs_c)
    return match.group(0) if match else None


def _find_matching_transcript_consequence(result: dict, transcript: str,
                                           debug: bool = False) -> Optional[dict]:
    """
    VEP can return consequences for MANY overlapping transcripts of the
    same gene for a single query -- confirmed live for P3H1, which
    returned 49 entries for one query, NINE of which shared the exact
    same coding notation (c.1838+2T>C) across different Ensembl
    transcript IDs. The `hgvsc` field on each entry uses Ensembl's OWN
    internal ENST accession, NEVER the originally-queried RefSeq
    accession (even though the top-level response echoes the RefSeq
    query back) -- confirmed live, so matching on hgvsc-startswith
    (the old approach) can never succeed when multiple transcripts are
    returned.

    This tries several PLAUSIBLE field names that Ensembl's confirmed,
    real `xref_refseq` VEP option might add to cross-reference each
    entry back to its RefSeq transcript ID(s) -- the exact JSON field
    name for this option could not be confirmed from documentation
    alone (unlike the "intron"/"exon" NUMBER/TOTAL format, which IS
    confirmed from Ensembl's own docs). If none of these match, this
    prints the FULL raw structure of every entry so a live run gives
    real ground truth to finalize this against, rather than guessing a
    fourth field name blind.
    """
    consequences = result.get("transcript_consequences", [])

    # Candidate field names that might hold RefSeq cross-reference IDs,
    # tried in order. Each is checked as either a single string or a list.
    candidate_fields = ("refseq_transcript_ids", "xref_refseq", "refseq", "source_refseq")

    for tc in consequences:
        for field in candidate_fields:
            value = tc.get(field)
            if value is None:
                continue
            values = value if isinstance(value, list) else [value]
            if any(str(v).startswith(transcript) for v in values):
                if debug:
                    print(f"[ensembl_hgvs_source] matched via {field!r}={value!r}")
                return tc

    # Fallback: the old hgvsc-prefix check, kept in case some future
    # response (or a gene with fewer overlapping transcripts) DOES echo
    # the RefSeq accession in hgvsc directly.
    for tc in consequences:
        hgvsc = tc.get("hgvsc", "")
        if hgvsc.startswith(transcript):
            return tc

    if debug:
        print(f"[ensembl_hgvs_source] could not match any of {len(consequences)} "
              f"transcript_consequences entries to {transcript!r} via RefSeq "
              f"cross-reference fields {candidate_fields!r} or hgvsc prefix. "
              f"Dumping the FULL first 3 raw entries so the real field name/shape "
              f"can be identified:")
        for i, tc in enumerate(consequences[:3]):
            print(f"[ensembl_hgvs_source]   entry #{i}: {tc!r}")
    return None


def get_splice_site_context(
    hgvs_c_with_transcript: str,
    genome_build: str = "GRCh38",
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
):
    """
    Determine which exon boundary a splice-region variant sits at.

    Only meaningful for intronic offset notation (c.NNNN+k / c.NNNN-k).
    Returns an empty SpliceSiteContext (all fields None) on any failure
    or missing data -- this is deliberately non-raising, mirroring how
    SpliceAI lookups degrade gracefully, since this is a supplementary
    descriptive sentence, not a fact the rest of the pipeline depends on.
    """
    from models import SpliceSiteContext  # local import to avoid a cycle

    transcript, hgvs_c = _parse_transcript_and_hgvs(hgvs_c_with_transcript)
    offset = _parse_hgvs_intronic_offset(hgvs_c)
    if offset is None:
        if debug:
            print(f"[ensembl_hgvs_source] {hgvs_c!r} has no intronic offset -- "
                  f"not a splice-region variant, skipping splice-context lookup")
        return SpliceSiteContext()

    base_url = _base_url_for_build(genome_build)
    http = session or requests.Session()
    url = f"{base_url}/vep/human/hgvs/{quote(hgvs_c_with_transcript, safe='')}"
    params = {"content-type": "application/json", "numbers": "1", "hgvs": "1", "xref_refseq": "1"}

    if debug:
        print(f"[ensembl_hgvs_source] GET {url} params={params} (splice-context lookup)")

    try:
        response = _get_with_retries(http, url, params, timeout_seconds, debug=debug)
    except requests.RequestException as exc:
        if debug:
            print(f"[ensembl_hgvs_source] splice-context request failed after retries: {exc}")
        return SpliceSiteContext()

    if debug:
        print(f"[ensembl_hgvs_source] splice-context status {response.status_code}")
        print(f"[ensembl_hgvs_source] splice-context response (first 800 chars): "
              f"{response.text[:800]}")

    if response.status_code != 200:
        return SpliceSiteContext()

    try:
        payload = response.json()
    except ValueError:
        return SpliceSiteContext()

    if not isinstance(payload, list) or not payload:
        return SpliceSiteContext()

    result = payload[0]
    tc = _find_matching_transcript_consequence(result, transcript, debug=debug)
    if tc is None:
        return SpliceSiteContext()

    intron_field = tc.get("intron")
    if not intron_field or "/" not in intron_field:
        if debug:
            print(f"[ensembl_hgvs_source] matched transcript_consequences entry has no "
                  f"usable 'intron' field: {intron_field!r} -- did the 'numbers' param "
                  f"actually turn this on? Full entry: {tc!r}")
        return SpliceSiteContext()

    try:
        intron_number = int(intron_field.split("/")[0])
    except ValueError:
        return SpliceSiteContext()

    is_donor = offset.startswith("+")
    if is_donor:
        # Donor site: intron immediately FOLLOWS this exon.
        site_type = "donor"
        nearest_exon = intron_number
    else:
        # Acceptor site: intron immediately PRECEDES the next exon.
        site_type = "acceptor"
        nearest_exon = intron_number + 1

    if debug:
        print(f"[ensembl_hgvs_source] splice context: offset={offset}, "
              f"intron={intron_number}, site_type={site_type}, nearest_exon={nearest_exon}")

    return SpliceSiteContext(offset=offset, site_type=site_type, nearest_exon=nearest_exon)
