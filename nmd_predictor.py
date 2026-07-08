"""
Nonsense-mediated decay (NMD) prediction for premature termination codons
(nonsense and frameshift variants).

RULE IMPLEMENTED (a SIMPLIFIED version of the standard Nagy & Maquat
1998 rule, per explicit request -- the full rule additionally requires
the PTC to be within the last ~50bp specifically of the second-to-last
exon, not just anywhere in it; that refinement is intentionally NOT
applied here):

    - PTC in the LAST exon                    -> escapes NMD
    - PTC in the second-to-last exon           -> escapes NMD
      (anywhere in it, not just its last ~50bp)
    - Everywhere else (any earlier exon)       -> triggers NMD

This is a real simplification, not just a wording change: the full rule
and this one only disagree for PTCs in the FIRST PART of the
second-to-last exon (full rule: triggers; this rule: escapes). They
agree everywhere else, including for a PTC in an early-to-middle exon
(e.g. exon 10 of 23) -- that still correctly triggers NMD under either
version.

WHY THIS WORKS IN cDNA (WHOLE-TRANSCRIPT) COORDINATES THROUGHOUT: exon
genomic lengths naturally sum to whole-transcript (cDNA) length,
including UTRs. Rather than trying to reconcile CDS-only boundaries with
genomic exon spans (which requires knowing exactly where the start/stop
codons fall within their exons), everything here is computed in cDNA
coordinates: the PTC's cDNA position comes directly from VEP, and the
final exon-exon junction's cDNA position comes from summing exon
lengths in transcript order. Both quantities being in the same
coordinate system is all that's needed for a correct relative-distance
comparison -- the absolute coordinate system doesn't matter. The
distance is still computed and reported (for reference/debugging) even
though it no longer gates the tier-2 decision.

CONFIRMED (from Ensembl's own official documentation):
  - GET /overlap/id/:id supports feature=exon, returning features that
    overlap a transcript, including its exons.

NOT YET LIVE-VERIFIED (flagging honestly, same as the earlier
xref_refseq situation): the exact JSON field names for each exon entry
(assumed candidates: "start"/"end" for genomic coordinates, "rank" for
transcript-order position -- "rank" is standard, well-established
Ensembl Core API terminology for a feature's position within its parent
transcript). This module checks several plausible field names and, if
none match, prints the full raw response so the real schema can be
confirmed from one live run rather than guessed further.

FRAMESHIFT PTC POSITION IS AN APPROXIMATION: for a frameshift variant,
VEP gives us the position where the frameshift BEGINS (cdna_start), not
where the new reading frame's stop codon actually falls. The PTC's cDNA
position is approximated as cdna_start + (stop_offset * 3), where
stop_offset is the number after "fs*"/"fsTer" in the HGVS p. notation
(e.g. 5 in "p.Leu1974ArgfsTer5"). This is the standard approximation
used for this kind of calculation (3 nucleotides per amino acid in the
shifted frame) but is still an approximation, not an exact genomic
lookup of the new stop codon -- flagging this so it isn't mistaken for
a directly-measured value.
"""
from __future__ import annotations
from typing import Optional
import re
import requests
from urllib.parse import quote

from ensembl_hgvs_source import (
    _base_url_for_build, _get_with_retries, _parse_transcript_and_hgvs,
    _find_matching_transcript_consequence, DEFAULT_TIMEOUT_SECONDS,
)
from models import NMDPrediction

NMD_BOUNDARY_NT = 50  # the standard "~50 bp rule" boundary

_FRAMESHIFT_STOP_OFFSET_PATTERN = re.compile(
    r"p\.\(?[A-Za-z]{3}\d+[A-Za-z]*fs(?:Ter|\*)(?P<stop_offset>\d+)"
)


def _parse_frameshift_stop_offset(hgvs_p: str) -> Optional[int]:
    """Extract the 'N' in 'fsTer N' / 'fs*N' from frameshift HGVS p.
    notation, e.g. 5 from 'p.Leu1974ArgfsTer5'. Returns None if hgvs_p
    doesn't match this pattern (e.g. no predicted stop position given)."""
    match = _FRAMESHIFT_STOP_OFFSET_PATTERN.search(hgvs_p or "")
    return int(match.group("stop_offset")) if match else None


def _fetch_ptc_position_data(
    hgvs_c_with_transcript: str,
    genome_build: str = "GRCh38",
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> Optional[dict]:
    """
    Query Ensembl VEP for the PTC-relevant position data: which exon the
    variant falls in, the total exon count, and the variant's own cDNA
    position. Reuses the exact same verified query parameters and
    RefSeq-cross-reference transcript-matching logic already proven for
    splice-site context lookups.

    Returns a dict with keys {exon_number, total_exons, cdna_start,
    ensembl_transcript_id} or None on any failure -- never guesses.
    """
    transcript, hgvs_c = _parse_transcript_and_hgvs(hgvs_c_with_transcript)
    base_url = _base_url_for_build(genome_build)
    http = session or requests.Session()
    url = f"{base_url}/vep/human/hgvs/{quote(hgvs_c_with_transcript, safe='')}"
    params = {"content-type": "application/json", "numbers": "1", "hgvs": "1", "xref_refseq": "1"}

    if debug:
        print(f"[nmd_predictor] GET {url} params={params} (PTC position lookup)")

    try:
        response = _get_with_retries(http, url, params, timeout_seconds, debug=debug)
    except requests.RequestException as exc:
        if debug:
            print(f"[nmd_predictor] PTC position request failed: {exc}")
        return None

    if debug:
        print(f"[nmd_predictor] PTC position status {response.status_code}")
        print(f"[nmd_predictor] PTC position response (first 800 chars): {response.text[:800]}")

    if response.status_code != 200:
        return None

    try:
        payload = response.json()
    except ValueError:
        return None
    if not isinstance(payload, list) or not payload:
        return None

    result = payload[0]
    tc = _find_matching_transcript_consequence(result, transcript, debug=debug)
    if tc is None:
        return None

    exon_field = tc.get("exon")
    cdna_start = tc.get("cdna_start")
    ensembl_transcript_id = tc.get("transcript_id")

    if not exon_field or "/" not in exon_field or cdna_start is None or not ensembl_transcript_id:
        if debug:
            print(f"[nmd_predictor] matched transcript_consequences entry is missing "
                  f"required fields (exon/cdna_start/transcript_id). Full entry: {tc!r}")
        return None

    try:
        exon_number, total_exons = (int(x) for x in exon_field.split("/"))
        cdna_start = int(cdna_start)
    except ValueError:
        return None

    if debug:
        print(f"[nmd_predictor] matched entry: ensembl_transcript_id={ensembl_transcript_id}, "
              f"exon={exon_number}/{total_exons}, cdna_start={cdna_start}, "
              f"refseq_transcript_ids={tc.get('refseq_transcript_ids')!r}")
        print(f"[nmd_predictor] IMPORTANT: if refseq_transcript_ids above lists more than "
              f"just {transcript!r}, this Ensembl transcript model is shared across "
              f"multiple RefSeq accessions and may not be an exact base-for-base match "
              f"to {transcript!r}'s own specific exon structure -- if the final NMD call "
              f"looks wrong, this is the first thing to suspect.")

    return {
        "exon_number": exon_number,
        "total_exons": total_exons,
        "cdna_start": cdna_start,
        "ensembl_transcript_id": ensembl_transcript_id,
    }


def _fetch_transcript_exon_lengths(
    ensembl_transcript_id: str,
    genome_build: str = "GRCh38",
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> Optional[list[int]]:
    """
    Fetch the genomic length of every exon in a transcript, ordered by
    transcript position (5' to 3' of the mature mRNA) -- exon[0] is the
    first exon, exon[-1] is the last.

    CONFIRMED endpoint (Ensembl's own overlap/id documentation):
        GET /overlap/id/{transcript_id}?feature=exon

    CONFIRMED REAL BUG FIX: overlap/id retrieves features that overlap
    the GENOMIC REGION spanned by the given ID -- it is NOT scoped to
    "only exons that belong to this transcript". For a gene with
    multiple transcript isoforms (or an overlapping neighboring gene),
    this can return exons from OTHER transcripts sharing that genomic
    space. Confirmed via a live case: querying BRCA1 (which has several
    annotated isoforms) without filtering produced a computed
    exon-junction distance of over 261,000 bp -- more than 30x the
    entire ~7.8kb BRCA1 mRNA length, an impossible value that could only
    come from summing exons across multiple different transcripts.

    Fix: filter the returned exon list to keep ONLY entries whose
    "Parent" field (confirmed as the standard field Ensembl's overlap
    endpoints use to identify a feature's owning parent -- seen directly
    on transcript-feature responses from this same endpoint family,
    e.g. {'Parent': 'ENSG00000234566', ...} for a transcript's parent
    gene; by the identical GFF3-style convention, an exon's Parent is
    its owning transcript) matching our exact target transcript ID
    BEFORE sorting or summing anything.

    NOT YET LIVE-VERIFIED: whether "Parent" is exactly the field name
    for exon-to-transcript ownership specifically (confirmed by analogy
    to a transcript-to-gene example, not a direct exon example), and
    whether it's a bare ID string or a version-qualified one (may need
    prefix-matching against the transcript ID rather than exact
    equality, to tolerate a versioned Parent like "ENST00000236040.9"
    against an unversioned query ID). Also tries "rank" (standard
    Ensembl terminology for transcript-order position) to sort directly
    within the filtered set; if absent, falls back to sorting by genomic
    start position and reversing for a minus-strand transcript. Prints
    the full raw response, and the filtered/unfiltered counts, so the
    real schema can be confirmed and any remaining discrepancy
    diagnosed from a live run.
    """
    base_url = _base_url_for_build(genome_build)
    http = session or requests.Session()
    url = f"{base_url}/overlap/id/{quote(ensembl_transcript_id, safe='')}"
    params = {"feature": "exon", "content-type": "application/json"}

    if debug:
        print(f"[nmd_predictor] GET {url} params={params} (fetching exon structure)")

    try:
        response = _get_with_retries(http, url, params, timeout_seconds, debug=debug)
    except requests.RequestException as exc:
        if debug:
            print(f"[nmd_predictor] exon structure request failed: {exc}")
        return None

    if debug:
        print(f"[nmd_predictor] exon structure status {response.status_code}")
        print(f"[nmd_predictor] exon structure response (first 1000 chars): {response.text[:1000]}")

    if response.status_code != 200:
        return None

    try:
        exons = response.json()
    except ValueError:
        return None
    if not isinstance(exons, list) or not exons:
        if debug:
            print(f"[nmd_predictor] unexpected exon list response shape: {exons!r}")
        return None

    # CRITICAL FILTER: keep only exons actually belonging to our target
    # transcript. Matches by prefix (not exact equality) to tolerate a
    # versioned Parent field (e.g. "ENST00000236040.9") against an
    # unversioned query ID, or vice versa.
    unfiltered_count = len(exons)
    own_exons = [e for e in exons if str(e.get("Parent", "")).startswith(ensembl_transcript_id)
                 or ensembl_transcript_id.startswith(str(e.get("Parent", "")) or "\0")]
    if debug:
        print(f"[nmd_predictor] {len(own_exons)}/{unfiltered_count} returned exons belong "
              f"to {ensembl_transcript_id} after Parent-field filtering")

    if not own_exons:
        if debug:
            print(f"[nmd_predictor] no exons matched {ensembl_transcript_id} via 'Parent' "
                  f"field -- cannot confirm which exons are ours. Raw Parent values seen: "
                  f"{[e.get('Parent') for e in exons]!r}. Refusing to guess by using the "
                  f"unfiltered (possibly multi-transcript) list.")
        return None

    exons = own_exons

    try:
        has_rank = all("rank" in e for e in exons)
        if has_rank:
            ordered = sorted(exons, key=lambda e: e["rank"])
        else:
            # Fallback: sort by genomic start, reverse for minus strand.
            if debug:
                print("[nmd_predictor] no 'rank' field on exon entries -- falling back "
                      "to genomic-coordinate sorting (less certain than 'rank')")
            strand = exons[0].get("strand", 1)
            ordered = sorted(exons, key=lambda e: e["start"])
            if strand in (-1, "-1"):
                ordered = list(reversed(ordered))

        lengths = [int(e["end"]) - int(e["start"]) + 1 for e in ordered]
    except (KeyError, ValueError, TypeError) as exc:
        if debug:
            print(f"[nmd_predictor] could not parse exon start/end/rank fields "
                  f"({exc}). Full raw exon list: {exons!r}")
        return None

    if any(length <= 0 for length in lengths):
        if debug:
            print(f"[nmd_predictor] computed a non-positive exon length -- "
                  f"something is wrong with field interpretation. Lengths: {lengths!r}")
        return None

    return lengths


def _compute_last_junction_cdna_position(exon_lengths: list[int]) -> Optional[int]:
    """
    Given transcript-ordered exon lengths, return the cDNA position of
    the LAST base of the second-to-last exon (i.e. the boundary between
    the final two exons -- the "NMD boundary" reference point). Returns
    None if there are fewer than 2 exons (single-exon transcripts have
    no exon-exon junction and are a documented special case: NMD does
    not apply at all to single-exon transcripts, since there's no
    splicing-deposited exon junction complex to trigger it).
    """
    if len(exon_lengths) < 2:
        return None
    return sum(exon_lengths[:-1])  # cumulative length up to (not including) the last exon


def predict_nmd(
    hgvs_c_with_transcript: str,
    variant_type: str,
    hgvs_p: Optional[str] = None,
    genome_build: str = "GRCh38",
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> NMDPrediction:
    """
    Predict whether a nonsense/frameshift variant's premature
    termination codon (PTC) is expected to trigger or escape
    nonsense-mediated decay.

    Only meaningful for variant_type in ('nonsense', 'frameshift').
    Returns an empty NMDPrediction (escapes_nmd=None) for any other
    variant type, or if any required data couldn't be resolved --
    degrades gracefully, never guesses.
    """
    if variant_type not in ("nonsense", "frameshift"):
        return NMDPrediction()

    ptc_data = _fetch_ptc_position_data(
        hgvs_c_with_transcript, genome_build=genome_build,
        timeout_seconds=timeout_seconds, session=session, debug=debug,
    )
    if ptc_data is None:
        return NMDPrediction()

    exon_number = ptc_data["exon_number"]
    total_exons = ptc_data["total_exons"]

    # Single-exon transcripts: NMD doesn't apply (no exon-exon junction
    # exists to deposit the exon junction complex that triggers it).
    if total_exons < 2:
        return NMDPrediction(
            escapes_nmd=True,
            reason="the transcript has only a single exon, so no exon-exon "
                   "junction exists to trigger NMD",
            ptc_exon_number=exon_number, total_exon_count=total_exons,
        )

    # Tier 1: PTC in the last exon -> escapes NMD, no further lookup needed.
    if exon_number == total_exons:
        return NMDPrediction(
            escapes_nmd=True, reason="the premature termination codon occurs in the last exon",
            ptc_exon_number=exon_number, total_exon_count=total_exons,
        )

    # Need the PTC's cDNA position (nonsense: the variant's own position;
    # frameshift: approximated downstream by stop_offset * 3) and the
    # exon structure to compute distance to the final junction.
    ptc_cdna_position = ptc_data["cdna_start"]
    if variant_type == "frameshift":
        stop_offset = _parse_frameshift_stop_offset(hgvs_p)
        if stop_offset is None:
            if debug:
                print("[nmd_predictor] frameshift variant has no parseable stop_offset "
                      "in hgvs_p -- cannot locate the PTC, skipping NMD prediction")
            return NMDPrediction()
        ptc_cdna_position += stop_offset * 3

    exon_lengths = _fetch_transcript_exon_lengths(
        ptc_data["ensembl_transcript_id"], genome_build=genome_build,
        timeout_seconds=timeout_seconds, session=session, debug=debug,
    )
    if exon_lengths is None:
        return NMDPrediction()

    junction_position = _compute_last_junction_cdna_position(exon_lengths)
    if junction_position is None:
        return NMDPrediction()

    distance = junction_position - ptc_cdna_position

    if debug:
        print(f"[nmd_predictor] PTC cDNA position={ptc_cdna_position}, "
              f"last junction cDNA position={junction_position}, distance={distance}, "
              f"PTC exon={exon_number}/{total_exons} (tier check: last exon={total_exons}, "
              f"second-to-last={total_exons - 1})")

    # Tier 2 (SIMPLIFIED, per explicit request): PTC anywhere in the
    # second-to-last exon escapes NMD, regardless of exact distance to
    # the final junction. This drops the standard rule's ~50bp
    # refinement (which only treats the LAST 50bp of the penultimate
    # exon as escaping) in favor of a simpler "last two exons" rule some
    # labs use for a quick first-pass call. distance_to_last_junction is
    # still computed and reported for reference, but no longer GATES
    # the tier-2 decision.
    #
    # IMPORTANT, flagged directly rather than left implicit: this
    # simplification does NOT mean "any exon near the end" escapes --
    # exon_number must be EXACTLY total_exons - 1. A variant in an
    # earlier exon (e.g. exon 10 of 23) still correctly triggers NMD
    # under this rule, same as under the full standard rule -- the two
    # only differ for PTCs specifically inside the second-to-last exon.
    if exon_number == total_exons - 1:
        return NMDPrediction(
            escapes_nmd=True,
            reason=f"the premature termination codon occurs in the second-to-last exon "
                   f"({distance} bp from the final exon-exon junction)",
            ptc_exon_number=exon_number, total_exon_count=total_exons,
            distance_to_last_junction=distance,
        )

    # Tier 3: everywhere else -> triggers NMD.
    return NMDPrediction(
        escapes_nmd=False,
        reason=f"the premature termination codon occurs in exon {exon_number} of "
               f"{total_exons}, neither the last nor second-to-last exon "
               f"({distance} bp upstream of the final exon-exon junction)",
        ptc_exon_number=exon_number, total_exon_count=total_exons,
        distance_to_last_junction=distance,
    )
