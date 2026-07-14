"""
ClinVar lookup via NCBI E-utilities (esearch + esummary + efetch).

We go straight to NCBI rather than relying on MyVariant.info's ClinVar
mirror for two of your six criteria specifically:

  1. "What is reported by other labs" -- this needs per-submission
     (per-SCV) classification data, which lives in ClinVar's full VCV
     record, not in the single aggregated "clinical_significance" field
     most aggregators expose.
  2. PubMed IDs -- ClinVar already curates literature citations per
     variant, so this doubles as your "PubMed IDs needed" criterion
     without a separate PubMed search.

Two NCBI calls per variant:
  - esearch (db=clinvar) to resolve genomic coordinates -> Variation ID
  - efetch (db=clinvar, rettype=vcv) to get the full record: aggregate
    classification, per-submitter classifications, and citation PMIDs

IMPORTANT: ClinVar's VCV XML schema has changed across NCBI releases, and
field paths below were written against the documented schema as of this
writing. If NCBI updates the schema, the XPath expressions in
`parse_vcv_xml` are the first thing to re-check against a live record.
"""
from __future__ import annotations
from typing import Optional
import re
import time
import xml.etree.ElementTree as ET
import requests

from models import ClinVarData, ClinVarSubmission, CoLocatedVariant

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TIMEOUT_SECONDS = 15# NCBI asks for no more than ~3 requests/sec without an API key.
REQUEST_DELAY_SECONDS = 0.4

def _normalize_hgvs_for_comparison(hgvs: str) -> Optional[tuple[str, str]]:
    """
    Extract (transcript, c._change) from either a plain HGVS string
    ("NM_000091.5:c.1408+2T>C") or a ClinVar VariationName
    ("NM_000091.5(COL4A3):c.1408+1G>C" or with a trailing protein
    consequence like "...(p.Gly2793Arg)"). Returns None if the string
    doesn't contain a recognizable transcript + c. description.

    CONFIRMED REAL BUG FIX: the deleted/duplicated sequence in del/dup
    HGVS notation is OPTIONAL -- "c.4065_4068del" and
    "c.4065_4068delTCAA" describe the exact same variant (ClinVar's own
    displayed VariationName omits it; some tools/labs spell it out).
    Without normalizing this away, a query with the sequence spelled out
    would never match ClinVar's own un-suffixed name for the identical
    variant -- confirmed live: querying "c.4065_4068delTCAA" against
    ClinVar's real "c.4065_4068del (p.Asn1355fs)" record for BRCA1 was
    rejected as a non-match, even though it's the correct record.

    The trailing-sequence strip only fires on a BARE "del"/"dup" ending
    (anchored to end-of-string) -- it deliberately does NOT touch
    "delins..." (a deletion-insertion), where the sequence after "ins"
    is NOT optional and describes what was actually inserted. The regex
    is naturally safe here: right after "del" in "delinsATG" comes "i",
    which isn't a valid base character, so the strip pattern simply
    doesn't match at that position.
    """
    if not hgvs:
        return None
    transcript_match = re.match(r"^([A-Za-z0-9_]+\.\d+)", hgvs.strip())
    change_match = re.search(r"(c\.[^\s(]+)", hgvs)
    if not transcript_match or not change_match:
        return None
    change = change_match.group(1).rstrip(".,;")
    change = re.sub(r"(del|dup)[ACGTacgt]+$", r"\1", change)
    return transcript_match.group(1), change


def _hgvs_matches(candidate_description: str, query_hgvs: str) -> bool:
    """True only if both the transcript accession AND the exact c. change
    match -- a shared transcript with a different position/change (the
    live bug this fixes: querying c.1408+2T>C matched a real ClinVar
    record for c.1408+1G>C) must be rejected, not accepted."""
    a = _normalize_hgvs_for_comparison(candidate_description)
    b = _normalize_hgvs_for_comparison(query_hgvs)
    if a is None or b is None:
        return False
    return a == b


def find_variation_id_by_hgvs(hgvs: str, debug: bool = False) -> Optional[str]:
    """
    Resolve an HGVS expression (e.g. NM_022356.4:c.1838+2T>C) directly to
    a ClinVar Variation ID.

    IMPORTANT -- verified against a real failure: NCBI's esearch does NOT
    do exact string matching on HGVS text. It splits the query into
    independent phrase terms (transcript, c. notation) and ANDs them,
    which means it can return a DIFFERENT, nearby variant that happens to
    share enough of those phrases. In a live run, searching for
    "NM_000091.5:c.1408+2T>C" matched a real ClinVar record whose actual
    VariationName was "NM_000091.5(COL4A3):c.1408+1G>C" -- a different
    intronic position (+1 vs +2) and a different substitution (G>C vs
    T>C). The earlier version of this function returned that first ID
    with no verification at all, which would have silently produced a
    complete-looking blurb (condition, PMIDs, classification, closing
    summary) describing the WRONG variant.

    This version fetches each candidate's VCV record and checks that its
    VariationName's HGVS description actually matches the query before
    accepting it -- the same stale-cache defense already used for
    VariantValidator responses (see hgvs_source_spare.py's
    _verify_match). Returns None (triggering the coordinate-based
    fallback in lookup_clinvar) if no candidate can be confirmed to
    match, rather than trusting an unverified text search.
    """
    resp = _ncbi_get(
        "esearch.fcgi",
        {"db": "clinvar", "term": hgvs, "retmode": "json"},
        debug=debug,
    )
    if resp is None:
        return None

    try:
        ids = resp.json()["esearchresult"]["idlist"]
    except (KeyError, ValueError):
        return None

    if debug:
        print(f"[clinvar_source] HGVS search: {hgvs}")
        print(f"[clinvar_source] candidate variation IDs (unverified): {ids}")

    for candidate_id in ids:
        root = fetch_vcv_xml(candidate_id, debug=debug)
        if root is None:
            continue
        vcv = root.find(".//VariationArchive")
        variation_name = vcv.get("VariationName") if vcv is not None else None
        if variation_name and _hgvs_matches(variation_name, hgvs):
            if debug:
                print(f"[clinvar_source] verified match: {candidate_id} "
                      f"(VariationName {variation_name!r} matches query {hgvs!r})")
            return candidate_id
        elif debug:
            print(f"[clinvar_source] REJECTED candidate {candidate_id}: "
                  f"VariationName {variation_name!r} does not match query {hgvs!r}")

    if debug and ids:
        print(f"[clinvar_source] no candidate could be verified as matching "
              f"{hgvs!r} -- falling back to coordinate-based search instead "
              f"of trusting an unverified text match.")
    return None


_REF_AA_POSITION_PATTERN = re.compile(r"p\.\(?([A-Za-z]{3})(\d+)")


def _extract_ref_aa_and_position(text: str) -> Optional[tuple[str, int]]:
    """
    Extract (reference_amino_acid, position) from any string containing
    HGVS protein notation -- works on a bare hgvs_p ('p.Asn1355Lysfs*10')
    or a full ClinVar VariationName containing a trailing parenthetical
    ('...(p.Asn1355fs)'), since this searches anywhere in the string
    rather than anchoring to the start.

    This works reliably across variant TYPES (missense, nonsense,
    frameshift all start with the reference amino acid + position, by
    HGVS definition) because the reference amino acid at a given codon
    is a fixed, definitional fact of the gene sequence -- any variant
    affecting that residue is described starting with the same prefix,
    regardless of what it changes to.
    """
    match = _REF_AA_POSITION_PATTERN.search(text or "")
    if not match:
        return None
    return match.group(1), int(match.group(2))


_BENIGN_ONLY_LABELS = {"benign", "likely benign"}


def _is_benign_only(classification: Optional[str]) -> bool:
    """
    True if a classification is PURELY within the benign spectrum --
    "Benign", "Likely benign", or a compound of just those two (e.g.
    "Benign/Likely benign"). Used to filter co-located variants that
    aren't usually clinically noteworthy to report alongside a variant
    being classified.

    A None/empty classification returns False (NOT considered
    benign-only) -- we don't want to silently drop a co-located variant
    just because we don't know its status; only a CONFIRMED benign-only
    call should be filtered out.

    Splits on "/" so compound labels are checked part-by-part (matching
    how classifications are combined elsewhere in this pipeline, e.g.
    "Pathogenic/Likely pathogenic") -- a compound like "Uncertain
    significance/Likely benign" is NOT benign-only (mixed), so it's kept.
    """
    if not classification:
        return False
    parts = [p.strip().lower() for p in classification.split("/") if p.strip()]
    return bool(parts) and all(p in _BENIGN_ONLY_LABELS for p in parts)


def find_co_located_variants(
    transcript: str,
    hgvs_p: str,
    exclude_variation_id: Optional[str] = None,
    max_results: int = 3,
    exclude_benign_only: bool = True,
    debug: bool = False,
) -> list[CoLocatedVariant]:
    """
    Find OTHER ClinVar variants affecting the SAME amino acid position as
    the one being reported (e.g. query is p.Asn1355Lysfs*10; this finds
    a separately-reported p.Asn1355Lys at the same residue) -- confirmed
    against a real reference report: "Another variant in the same
    position, p.Arg5105Gln, has conflicting classifications of
    pathogenicity by other clinical laboratories (ClinVar ID: 69441)."

    APPROACH: search ClinVar's text index for {transcript} + the
    reference-amino-acid+position prefix (e.g. "NM_007294.4 p.Asn1355"),
    the same free-text search mechanism already confirmed to work (and
    to require verification -- see find_variation_id_by_hgvs's docstring
    for the exact live case that proved NCBI's esearch does phrase-level
    matching, not exact matching). Every candidate is fetched and its
    OWN position is re-extracted and compared before being trusted --
    this is NOT an exact-match search (that's the whole point -- we WANT
    different changes at the same position), but the POSITION NUMBER
    itself must match exactly, and the transcript must be the same one
    being reported against.

    SCOPE, stated explicitly: only searches within this exact transcript
    accession -- deliberately does not attempt to also find records
    submitted against a different transcript version of the same gene,
    consistent with this pipeline's "always verify against the exact
    requested transcript" approach used everywhere else.

    Args:
        transcript: e.g. "NM_007294.4" -- the transcript being reported against.
        hgvs_p: the QUERY variant's own protein notation, used to derive
            the reference-amino-acid+position search term.
        exclude_variation_id: if the query variant itself has a known
            ClinVar ID, pass it here so it doesn't get listed as
            "another" variant (it would otherwise match itself, since
            its own position obviously matches).
        max_results: cap on how many co-located variants to return, to
            avoid an overly long blurb if a position happens to have
            many reported variants (same capping pattern used for PMIDs
            elsewhere in this pipeline). Applied AFTER the benign filter
            below, so it caps the clinically-relevant set, not a mix
            that includes entries about to be discarded anyway.
        exclude_benign_only: if True (the default), skip any co-located
            variant whose classification is PURELY within the benign
            spectrum (Benign, Likely benign, or a compound of just
            those two, e.g. "Benign/Likely benign") -- these aren't
            usually clinically noteworthy to report alongside a variant
            being classified. A variant with NO known classification is
            NOT excluded by this (we don't want to silently drop a
            co-located variant just because we don't know its status --
            only a CONFIRMED benign-only call is filtered out).

    Returns:
        A list of CoLocatedVariant, ordered as returned by ClinVar's
        search (not otherwise sorted/ranked) -- empty if hgvs_p has no
        parseable position, the search fails, or nothing else was found
        at that position.
    """
    query_ref_aa_pos = _extract_ref_aa_and_position(hgvs_p)
    if query_ref_aa_pos is None:
        if debug:
            print(f"[clinvar_source] could not extract a reference-amino-acid+position "
                  f"prefix from hgvs_p={hgvs_p!r} -- skipping co-located variant search")
        return []
    query_ref_aa, query_position = query_ref_aa_pos

    term = f"{transcript} p.{query_ref_aa}{query_position}"
    resp = _ncbi_get("esearch.fcgi", {"db": "clinvar", "term": term, "retmode": "json"}, debug=debug)
    if resp is None:
        return []

    try:
        ids = resp.json()["esearchresult"]["idlist"]
    except (KeyError, ValueError):
        return []

    if debug:
        print(f"[clinvar_source] co-located variant search: {term!r}")
        print(f"[clinvar_source] candidate variation IDs (unverified): {ids}")

    results: list[CoLocatedVariant] = []
    for candidate_id in ids:
        if candidate_id == exclude_variation_id:
            if debug:
                print(f"[clinvar_source] skipping {candidate_id}: this is the query "
                      f"variant's own ClinVar record, not 'another' variant")
            continue

        root = fetch_vcv_xml(candidate_id, debug=debug)
        if root is None:
            continue
        vcv = root.find(".//VariationArchive")
        variation_name = vcv.get("VariationName") if vcv is not None else None
        candidate_ref_aa_pos = _extract_ref_aa_and_position(variation_name or "")

        if candidate_ref_aa_pos is None or candidate_ref_aa_pos[1] != query_position:
            if debug:
                print(f"[clinvar_source] REJECTED co-located candidate {candidate_id}: "
                      f"VariationName {variation_name!r} does not confirm position "
                      f"{query_position} -- likely an unrelated text match")
            continue

        # Extract just the protein-change portion for display, e.g.
        # "p.Asn1355Lys" from "...del (p.Asn1355Lys)".
        protein_match = re.search(r"p\.\(?[A-Za-z0-9*=]+\)?", variation_name)
        protein_change = protein_match.group(0).strip("()") if protein_match else variation_name

        # CONFIRMED BUG FIX: also extract the c. notation (same pattern
        # already used in _normalize_hgvs_for_comparison), needed to
        # disambiguate co-located variants whose protein notation
        # collapses to an identical string -- most notably synonymous
        # variants, where c.8379A>C, c.8379A>T, and c.8379A>G all
        # display as the same "p.Gly2793=" despite being three distinct
        # variants. Confirmed against a real live case (BRCA2
        # p.Gly2793Arg's co-located search returning three "p.Gly2793="
        # entries that were indistinguishable in the rendered sentence).
        hgvs_c_match = re.search(r"(c\.[^\s(]+)", variation_name)
        hgvs_c = hgvs_c_match.group(1).rstrip(".,;") if hgvs_c_match else None

        classified = parse_vcv_xml(root)
        if debug:
            print(f"[clinvar_source] confirmed co-located variant {candidate_id}: "
                  f"{protein_change!r} ({hgvs_c!r}), classification={classified.aggregate_classification!r}")

        if exclude_benign_only and _is_benign_only(classified.aggregate_classification):
            if debug:
                print(f"[clinvar_source] EXCLUDED co-located variant {candidate_id}: "
                      f"classification {classified.aggregate_classification!r} is purely "
                      f"benign/likely benign -- not reported (exclude_benign_only=True)")
            continue

        results.append(CoLocatedVariant(
            protein_change=protein_change,
            hgvs_c=hgvs_c,
            classification=classified.aggregate_classification,
            clinvar_id=candidate_id,
        ))
        if len(results) >= max_results:
            break

    return results


def _ncbi_get(endpoint: str, params: dict, debug: bool = False) -> Optional[requests.Response]:
    try:
        resp = requests.get(f"{EUTILS_BASE}/{endpoint}", params=params, timeout=TIMEOUT_SECONDS)
        time.sleep(REQUEST_DELAY_SECONDS)
        if debug:
            print(f"[clinvar_source] GET {endpoint} -> status {resp.status_code}, params={params}")
            print(f"[clinvar_source] response (first 3000 chars): {resp.text[:3000]}")
        if resp.status_code != 200:
            return None
        return resp
    except requests.RequestException as e:
        if debug:
            print(f"[clinvar_source] request failed: {e}")
        return None


def _to_vcv_accession(variation_id: str) -> str:
    """
    NCBI's efetch with rettype=vcv requires the VCV accession format
    (e.g. 'VCV000052569'), NOT the bare numeric Variation ID ('52569') --
    confirmed against NCBI's own documented example
    (https://ncbiinsights.ncbi.nlm.nih.gov/2019/10/11/visit-the-new-clinvar-for-easier-variant-interpretation/):
    'efetch.fcgi?db=clinvar&rettype=vcv&id=VCV000007105'

    Passing the bare number (as the original version of this module did)
    doesn't error -- it just silently returns an empty <ClinVarResult-Set>,
    which is exactly what broke this the first time.
    """
    digits = "".join(c for c in str(variation_id) if c.isdigit())
    return f"VCV{int(digits):09d}"


def fetch_vcv_by_variation_id(variation_id: str, debug: bool = False) -> Optional[ET.Element]:
    """
    Skips the coordinate-based esearch step entirely. Use this to sanity
    check the parsing logic against a Variation ID you already know (e.g.
    from a GeneBe/VarSeq export, or a ClinVar webpage), independent of
    whether your genomic coordinates are correct.
    """
    return fetch_vcv_xml(variation_id, debug=debug)


def find_candidate_variation_ids(chrom: str, pos: int, genome_build: str = "GRCh38",
                                  debug: bool = False) -> list[str]:
    """
    Resolve chr:pos to candidate ClinVar Variation IDs via esearch.

    Verified query syntax (confirmed against NCBI's own help docs and a
    working Bio.Entrez example, since the original [chromosome]/[Reference]/
    [Allele] field tags this module used were not real ClinVar search
    fields and were silently ignored by NCBI -- see the chrom/chrpos
    examples at https://www.ncbi.nlm.nih.gov/clinvar/docs/help/):
        17[chr] AND 43000000:44000000[chrpos37]

    Position is searched as a RANGE even for a single base. A given
    position can have more than one variant recorded (different ref>alt),
    so this returns ALL matches -- disambiguation by exact ref/alt happens
    in lookup_clinvar(), not here.
    """
    chrom_clean = chrom.replace("chr", "")
    build_tag = "38" if genome_build.upper() in ("GRCH38", "HG38") else "37"
    term = f"{chrom_clean}[chr] AND {pos}:{pos}[chrpos{build_tag}]"
    resp = _ncbi_get("esearch.fcgi", {"db": "clinvar", "term": term, "retmode": "json"}, debug=debug)
    if resp is None:
        return []
    try:
        ids = resp.json().get("esearchresult", {}).get("idlist", [])
        if debug:
            print(f"[clinvar_source] esearch term: {term}")
            print(f"[clinvar_source] matched variation IDs: {ids}")
        return ids
    except (ValueError, KeyError):
        return []


def _vcv_matches_allele(root: ET.Element, ref: str, alt: str, debug: bool = False) -> bool:
    """
    Best-effort check of whether a fetched VCV record's reported alleles
    match the ref/alt you're looking for, used to disambiguate when a
    single position has multiple ClinVar records.

    NOTE: this checks a couple of plausible attribute locations
    (ReferenceAlleleVCF/AlternateAlleleVCF on SequenceLocation) based on
    the documented VCV schema, but has not been confirmed against a live
    multi-variant-at-one-position record. If you hit a position with
    several variants and the wrong one gets selected, print the raw XML
    (debug=True) and check the actual attribute names/values, then adjust
    this function accordingly.
    """
    for loc in root.findall(".//SequenceLocation"):
        loc_ref = loc.get("referenceAlleleVCF") or loc.get("ReferenceAlleleVCF")
        loc_alt = loc.get("alternateAlleleVCF") or loc.get("AlternateAlleleVCF")
        if loc_ref and loc_alt:
            if debug:
                print(f"[clinvar_source] candidate alleles found in record: {loc_ref}>{loc_alt}")
            return loc_ref.upper() == ref.upper() and loc_alt.upper() == alt.upper()
    # If we can't find allele info to check, don't silently claim a match --
    # let the caller know this couldn't be verified.
    return False


def find_variation_id(chrom: str, pos: int, ref: str, alt: str, genome_build: str = "GRCh38",
                       debug: bool = False) -> Optional[str]:
    """
    Resolve chr:pos:ref:alt to a single ClinVar Variation ID, verifying
    the allele match against EVERY candidate -- including when there's
    only one.

    CONFIRMED BUG FIX: the previous version returned a lone candidate
    immediately, with NO allele verification, on the assumption that
    "only one candidate at this position" meant "this must be the right
    one." That's false: a position can have exactly one ClinVar record,
    and it can still be for a DIFFERENT substitution than the one being
    queried. Confirmed against a live case: querying TTN
    2:178747087 G>A (a nonsense variant, c.15313C>T/p.Arg5105*) returned
    ClinVar's ONLY record at that position, VCV000497554 -- but that
    record is for G>T (a different, synonymous variant, c.15313C>A/
    p.Arg5105=). The single-candidate shortcut skipped allele
    verification entirely and returned data for the wrong variant.

    Now returns None (i.e. "no verified match" -- callers treat this the
    same as "not in ClinVar") if the only candidate's alleles don't
    match, rather than silently attaching an unrelated variant's
    classification, PMIDs, and condition to the report.
    """
    candidates = find_candidate_variation_ids(chrom, pos, genome_build, debug=debug)
    if not candidates:
        return None

    for vid in candidates:
        root = fetch_vcv_xml(vid, debug=debug)
        if root is not None and _vcv_matches_allele(root, ref, alt, debug=debug):
            return vid

    if debug:
        print(f"[clinvar_source] {len(candidates)} candidate(s) found at this position "
              f"({candidates}), but NONE could be confirmed to match {ref}>{alt} -- "
              f"treating as 'no verified ClinVar match' rather than guessing. This is "
              f"correct behavior when the position has a DIFFERENT variant recorded "
              f"(e.g. a different substitution at the same coordinate) but not the "
              f"exact one being queried.")
    return None


def fetch_vcv_xml(variation_id: str, debug: bool = False) -> Optional[ET.Element]:
    vcv_accession = _to_vcv_accession(variation_id)
    resp = _ncbi_get("efetch.fcgi", {"db": "clinvar", "id": vcv_accession, "rettype": "vcv"}, debug=debug)
    if resp is None:
        return None
    try:
        return ET.fromstring(resp.text)
    except ET.ParseError as e:
        if debug:
            print(f"[clinvar_source] XML parse failed: {e}")
        return None


# ClinVar's generic placeholder condition names -- when no real disease
# was submitted. We skip these so the opening line doesn't read
# "reported in individuals with not provided".
_GENERIC_CONDITIONS = {
    "not provided", "not specified", "see cases", "none provided",
    "not applicable", "conditions", "allhighlypenetrant",
}


def _first_match_in(element: ET.Element, *xpaths: str) -> Optional[ET.Element]:
    """
    Return the first element matching any of the given XPaths, searched
    relative to `element`. Shared by both the aggregate classification
    lookup (relative to the document root) and the per-submission lookup
    (relative to each ClinicalAssertion), since both need to tolerate the
    same 2024+ VCV schema variations.

    IMPORTANT: uses explicit `is not None` checks, never `a or b` --  an
    ElementTree element with no children is falsy even if it has real
    text or attributes, so an or-chain silently skips a genuine match.
    """
    for xp in xpaths:
        el = element.find(xp)
        if el is not None:
            return el
    return None


def _extract_condition(root: ET.Element) -> Optional[str]:
    """
    Return a human-readable associated condition name, or None if only
    generic placeholders are present. Conditions live under TraitSet/Trait,
    with names in Name/ElementValue (VCV schema). Prefers a name marked
    Type='Preferred'; otherwise takes the first non-generic name found.
    """
    names: list[str] = []
    preferred: list[str] = []
    for trait in root.findall(".//TraitSet/Trait"):
        for name_el in trait.findall("./Name"):
            value_el = name_el.find("./ElementValue")
            if value_el is None or not value_el.text:
                continue
            text = value_el.text.strip()
            if text.lower() in _GENERIC_CONDITIONS:
                continue
            if value_el.get("Type") == "Preferred":
                preferred.append(text)
            else:
                names.append(text)
    if preferred:
        return preferred[0]
    if names:
        return names[0]
    return None


def parse_vcv_xml(root: ET.Element) -> ClinVarData:
    """
    Pulls: aggregate classification + review status, per-submitter
    classifications, and citation PMIDs from a VCV record.
    """
    data = ClinVarData()

    vcv = root.find(".//VariationArchive")
    if vcv is not None:
        data.variation_id = vcv.get("VariationID") or vcv.get("Accession")

    # Aggregate (overall) germline classification.
    # The 2024+ VCV XML redesign (which added somatic/oncogenicity
    # classifications) can nest the germline classification a couple of
    # different ways depending on record vintage. We check the known
    # layouts in order so the aggregate label resolves on both old and
    # current-format records. efetch returns the new format by default
    # (old format required old_xml=T, which we don't send), but records
    # vary, so we stay tolerant of both.
    agg = _first_match_in(
        root,
        ".//Classifications/GermlineClassification/Description",
        ".//GermlineClassification/Description",
        ".//ClinicalSignificance/Description",
    )
    if agg is not None:
        data.aggregate_classification = agg.text

    review = _first_match_in(
        root,
        ".//Classifications/GermlineClassification/ReviewStatus",
        ".//GermlineClassification/ReviewStatus",
        ".//ClinicalSignificance/ReviewStatus",
    )
    if review is not None:
        data.review_status = review.text

    # Associated condition/disease name, for the opening line. Conditions
    # live in TraitSet/Trait with a Name/ElementValue. ClinVar uses generic
    # placeholders ("not provided", "not specified") when no real disease
    # was submitted -- we skip those so the opening line doesn't say
    # "reported in individuals with not provided". Prefer a Preferred-type
    # name when available; fall back to the first real name otherwise.
    condition = _extract_condition(root)
    if condition:
        data.condition = condition

    # Per-submission ("other labs") classifications.
    # FIXED BUG: this used to assume GermlineClassification was always a
    # leaf text node, the same wrong assumption already corrected for the
    # AGGREGATE classification above. For the 2024+ VCV schema,
    # GermlineClassification can be a WRAPPER containing a <Description>
    # child -- reading .text directly on the wrapper returns only the
    # whitespace/newline before that child tag, which is a non-empty
    # (truthy) string, so a submission got appended with a blank
    # classification instead of being correctly skipped or resolved.
    # Confirmed against a live case: a single-submission ClinVar record
    # (VCV004340382) produced "This variant has been classified as  by
    # other clinical laboratories" -- a blank label -- via exactly this
    # path.
    for assertion in root.findall(".//ClinicalAssertion"):
        submitter_el = assertion.find(".//ClinVarAccession")
        submitter = submitter_el.get("SubmitterName") if submitter_el is not None else None

        classif_el = _first_match_in(
            assertion,
            ".//Classification/GermlineClassification/Description",
            ".//Classification/GermlineClassification",
            ".//Classification",
        )
        review_el = assertion.find(".//Classification/ReviewStatus")

        # Guard against whitespace-only text from a wrapper element with
        # no Description child matched -- treat that the same as "no
        # classification found" rather than a truthy-but-blank value.
        classification_text = classif_el.text.strip() if classif_el is not None and classif_el.text else None
        review_text = review_el.text if review_el is not None else None

        if submitter and classification_text:
            data.submissions.append(ClinVarSubmission(
                submitter=submitter,
                classification=classification_text,
                review_status=review_text,
            ))

    return data


def extract_pmids(root: ET.Element) -> list[str]:
    pmids = set()
    for citation in root.findall(".//Citation"):
        id_el = citation.find("./ID[@Source='PubMed']")
        if id_el is not None and id_el.text:
            pmids.add(id_el.text.strip())
    return sorted(pmids)


def lookup_clinvar(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    genome_build: str = "GRCh38",
    hgvs: str | None = None,
    debug: bool = False,
) -> tuple[ClinVarData, list[str]]:
    """Returns (ClinVarData, pubmed_ids)."""

    variation_id = None

    # First try HGVS
    if hgvs:
        variation_id = find_variation_id_by_hgvs(hgvs, debug=debug)

    # Fall back to coordinate lookup
    if variation_id is None:
        variation_id = find_variation_id(
            chrom, pos, ref, alt,
            genome_build=genome_build,
            debug=debug,
        )

    if variation_id is None:
        return ClinVarData(variation_id=None, status="not_found"), []

    root = fetch_vcv_xml(variation_id, debug=debug)
    if root is None:
        return ClinVarData(variation_id=variation_id, status="not_found"), []

    clinvar_data = parse_vcv_xml(root)

    if not clinvar_data.variation_id:
        clinvar_data.variation_id = variation_id

    if debug:
        print(f"[clinvar_source] parsed: aggregate_classification="
              f"{clinvar_data.aggregate_classification!r}, "
              f"{len(clinvar_data.submissions)} submission(s): "
              f"{[(s.submitter, s.classification) for s in clinvar_data.submissions]!r}")
        if not clinvar_data.aggregate_classification and not clinvar_data.submissions:
            # Nothing was extracted at all -- dump the raw structure of
            # every ClinicalAssertion found, so the actual tag/nesting
            # this record uses is visible instead of guessing again.
            assertions = root.findall(".//ClinicalAssertion")
            print(f"[clinvar_source] classification extraction found NOTHING for this "
                  f"record -- dumping raw structure of {len(assertions)} "
                  f"ClinicalAssertion element(s) found in the XML:")
            for i, assertion in enumerate(assertions):
                print(f"[clinvar_source]   assertion #{i}: "
                      f"{ET.tostring(assertion, encoding='unicode')[:2000]}")

    pmids = extract_pmids(root)

    return clinvar_data, pmids


def lookup_clinvar_by_variation_id(variation_id: str, debug: bool = False) -> tuple[ClinVarData, list[str]]:
    """
    Diagnostic helper: fetch and parse a known ClinVar Variation ID directly,
    bypassing coordinate resolution. Use this when you already know the ID
    (e.g. ClinVar ID 52569) and want to confirm the network/parsing layer
    works, independent of whether your genomic coordinates are correct.
    """
    root = fetch_vcv_by_variation_id(variation_id, debug=debug)
    if root is None:
        return ClinVarData(variation_id=variation_id, status="not_found"), []
    clinvar_data = parse_vcv_xml(root)
    if not clinvar_data.variation_id:
        clinvar_data.variation_id = variation_id
    return clinvar_data, extract_pmids(root)
