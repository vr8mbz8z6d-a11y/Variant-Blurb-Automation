"""
hgvs_source.py

Resolves transcript-level HGVS coding variant descriptions (e.g.
"NM_000162.5:c.1021G>A") into strand-corrected, VCF-style genomic
coordinates (chrom, pos, ref, alt) by calling the public VariantValidator
REST API.

Why this exists
----------------
Copying c. notation letters directly into VCF ref/alt only happens to be
correct for plus-strand genes. For minus-strand genes (e.g. GCK), the
c. notation is written relative to the transcript, not the genomic
plus strand, so ref/alt must be complemented. VariantValidator performs
this conversion correctly (it understands strand, splice-site offsets,
indel normalisation, and multi-exon spanning variants), so we wrap it
rather than reimplementing HGVS parsing ourselves.

Usage
-----
    from hgvs_source import resolve_coordinates, HGVSResolutionError

    result = resolve_coordinates("NM_000162.5:c.1021G>A", genome_build="GRCh38")
    # result.chrom, result.pos, result.ref, result.alt

    # Or, given only a gene symbol + c. notation (no transcript accession):
    from hgvs_source import resolve_coordinates_from_gene

    result = resolve_coordinates_from_gene("GCK", "c.1021G>A", genome_build="GRCh38")
    # result.transcript tells you which transcript was actually used

Design notes
------------
- Transcript selection: resolve_coordinates() requires the caller to
  supply the transcript accession as part of the hgvs_c string (e.g.
  "NM_000162.5:c.1021G>A", not just "c.1021G>A"), because which
  transcript to use is a clinically meaningful choice (different
  transcripts can number codons differently).
  resolve_coordinates_from_gene() offers a convenience path that picks
  the transcript automatically -- specifically the MANE Select
  transcript, which is the standardized single "primary" RefSeq
  transcript NCBI/EMBL-EBI designate per gene for exactly this purpose.
  It is NOT a silent guess among many equally-valid options: MANE Select
  is a defined, single answer per gene, and the resolved transcript is
  always returned on the result (`.transcript`) so it can be surfaced
  to a reviewer. Still prefer resolve_coordinates() with an explicit
  accession when you have one (e.g. from a VarSeq export), since it
  removes even this one API round-trip and any doubt about which
  transcript was used.
- Stale-cache guard: there is a known intermittent bug where the public
  VariantValidator REST API can return a cached result for a different
  query. To guard against this, resolve_coordinates() verifies that the
  returned variant description actually corresponds to the requested
  gene/transcript/HGVS notation before accepting the result, and raises
  HGVSResolutionError if it doesn't match.
- Fail closed, not silently: any failure (network error, non-2xx
  response, validation mismatch, unexpected payload shape) raises
  HGVSResolutionError. Callers with independently-known coordinates are
  expected to catch this and fall back to them, per the pipeline's
  existing behaviour. resolve_coordinates_from_gene() has no such
  fallback available -- gene + c. notation IS the only input -- so
  callers should let HGVSResolutionError propagate and surface it,
  rather than silently producing an empty result.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Optional
from urllib.parse import quote

import requests

DEFAULT_BASE_URL = "https://rest.variantvalidator.org"
DEFAULT_TIMEOUT_SECONDS = 30  # bumped from 15 after repeated live timeouts on
                              # splice-site queries; VariantValidator does real
                              # computation server-side, not a simple lookup, so
                              # 15s was cutting it close on some inputs. This may
                              # not fully fix transient server-side slowness, but
                              # it removes an easy, avoidable cause of failure.

# VariantValidator's API only recognises "GRCh37"/"GRCh38" as genome build
# values -- passing "hg38" straight through (a natural thing to do, since
# the rest of this pipeline uses "hg38"/"hg19" everywhere) doesn't error,
# it just silently returns a result with an empty primary_assembly_loci
# block, which is what produced "Available builds: []". We accept the
# common aliases here and normalize to what the API actually expects,
# once, so every caller (URL construction AND response parsing) agrees.
_GENOME_BUILD_ALIASES = {
    "hg38": "GRCh38",
    "grch38": "GRCh38",
    "38": "GRCh38",
    "hg19": "GRCh37",
    "grch37": "GRCh37",
    "37": "GRCh37",
}


def _normalize_genome_build(genome_build: str) -> str:
    """Map hg38/hg19/GRCh38/GRCh37 (any casing) to the canonical
    "GRCh38"/"GRCh37" string VariantValidator's API expects.
    """
    key = (genome_build or "").strip().lower()
    canonical = _GENOME_BUILD_ALIASES.get(key)
    if canonical is None:
        raise HGVSResolutionError(
            f"Unrecognised genome build {genome_build!r}. Expected one of: "
            f"hg19, hg38, GRCh37, GRCh38 (any casing)."
        )
    return canonical


_HGVS_WITH_TRANSCRIPT_RE = re.compile(
    r"^(?P<transcript>[A-Z]{1,3}_\d+\.\d+):(?P<hgvs_c>c\..+)$"
)


class HGVSResolutionError(Exception):
    """Raised when an HGVS description cannot be confidently resolved
    to genomic coordinates. Callers should catch this and fall back to
    manually-provided coordinates rather than trusting a partial or
    unverified result.
    """


@dataclass(frozen=True)
class ResolvedCoordinates:
    """Strand-corrected, VCF-style genomic coordinates for a variant."""

    chrom: str
    pos: int
    ref: str
    alt: str
    genome_build: str
    transcript: str
    hgvs_c: str
    # The full variant description string VariantValidator echoed back,
    # kept for audit/debug purposes.
    validated_description: str
    strand: Optional[int] = None


def _parse_transcript_and_hgvs(hgvs_c_with_transcript: str) -> tuple[str, str]:
    """Split "NM_000162.5:c.1021G>A" into ("NM_000162.5", "c.1021G>A").

    Raises HGVSResolutionError if the input doesn't include a transcript
    accession -- this module deliberately does not guess one.
    """
    match = _HGVS_WITH_TRANSCRIPT_RE.match(hgvs_c_with_transcript.strip())
    if not match:
        raise HGVSResolutionError(
            "Expected a transcript-qualified HGVS coding description like "
            "'NM_000162.5:c.1021G>A', got: "
            f"{hgvs_c_with_transcript!r}. This module does not auto-select "
            "a transcript -- pass the accession explicitly."
        )
    return match.group("transcript"), match.group("hgvs_c")


def _find_matching_result(payload: dict[str, Any], transcript: str) -> dict[str, Any]:
    """VariantValidator's response is keyed by variant description
    strings, one of which should correspond to our requested transcript.
    Return that entry, raising if none matches.
    """
    if not isinstance(payload, dict):
        raise HGVSResolutionError(
            f"Unexpected response shape from VariantValidator (not a dict): {type(payload)}"
        )

    candidate_keys = [k for k in payload.keys() if k not in ("metadata", "flag")]

    for key in candidate_keys:
        entry = payload.get(key)
        if isinstance(entry, dict) and transcript in key:
            return entry

    # Fall back: if there's exactly one real candidate, use it, but still
    # verify its contents below via _verify_match.
    if len(candidate_keys) == 1:
        entry = payload.get(candidate_keys[0])
        if isinstance(entry, dict):
            return entry

    raise HGVSResolutionError(
        f"Could not find a result for transcript {transcript!r} in "
        f"VariantValidator response. Keys returned: {candidate_keys!r}"
    )


def _verify_match(entry: dict[str, Any], transcript: str, hgvs_c: str) -> str:
    """Guard against the known stale-cache bug: confirm the returned
    entry actually corresponds to what we asked for before trusting its
    coordinates. Returns the validated description string used for the
    check (kept for audit purposes).

    Raises HGVSResolutionError if the entry doesn't correspond to the
    request.
    """
    description = (
        entry.get("hgvs_transcript_variant")
        or entry.get("submitted_variant")
        or ""
    )

    if not description:
        raise HGVSResolutionError(
            "VariantValidator response entry is missing a variant "
            "description to verify against the request; refusing to "
            "trust its coordinates (possible stale-cache response)."
        )

    if transcript not in description:
        raise HGVSResolutionError(
            f"Stale/mismatched VariantValidator response: expected transcript "
            f"{transcript!r} in returned description, got {description!r}. "
            "This matches a known intermittent caching bug -- refusing to "
            "trust these coordinates."
        )

    # The c. change itself (e.g. "c.1021G>A") should also appear, ignoring
    # whitespace differences.
    normalized_expected = hgvs_c.replace(" ", "")
    normalized_actual = description.replace(" ", "")
    if normalized_expected not in normalized_actual:
        raise HGVSResolutionError(
            f"Stale/mismatched VariantValidator response: expected HGVS change "
            f"{hgvs_c!r} in returned description, got {description!r}. "
            "This matches a known intermittent caching bug -- refusing to "
            "trust these coordinates."
        )

    return description


def _extract_vcf_fields(entry: dict[str, Any], genome_build: str, debug: bool = False) -> tuple[str, int, str, str]:
    """Pull chr/pos/ref/alt out of a VariantValidator result entry.

    VariantValidator nests the VCF block under the genome build key,
    e.g. entry["primary_assembly_loci"]["grch38"]["vcf"].
    """
    build_key = genome_build.lower()
    try:
        loci = entry["primary_assembly_loci"]
    except KeyError as exc:
        if debug:
            print(f"[hgvs_source] entry has no 'primary_assembly_loci' key at all. "
                  f"Entry keys: {list(entry.keys())!r}")
            print(f"[hgvs_source] validation_warnings: {entry.get('validation_warnings')!r}")
        raise HGVSResolutionError(
            "VariantValidator response entry has no 'primary_assembly_loci' field"
        ) from exc

    build_block = loci.get(build_key)
    if build_block is None:
        if debug:
            print(f"[hgvs_source] 'primary_assembly_loci' is present but has no "
                  f"{build_key!r} key. Full loci dict: {loci!r}")
            print(f"[hgvs_source] validation_warnings on this entry: "
                  f"{entry.get('validation_warnings')!r}")
            print(f"[hgvs_source] This most often means VariantValidator validated the "
                  f"HGVS description itself but could not project it onto primary "
                  f"GRCh38/GRCh37 assembly coordinates for this transcript -- a known "
                  f"issue category for genes with complex genome-transcript alignment "
                  f"(TTN is a well-known example due to its size and exon complexity). "
                  f"Check 'validation_warnings' above for the specific reason.")
        raise HGVSResolutionError(
            f"VariantValidator response has no results for genome build "
            f"{genome_build!r}. Available builds: {list(loci.keys())!r}"
        )

    vcf = build_block.get("vcf")
    if not vcf:
        if debug:
            print(f"[hgvs_source] build block for {build_key!r} present but has no "
                  f"'vcf' key. Full block: {build_block!r}")
        raise HGVSResolutionError(
            f"VariantValidator response for build {genome_build!r} has no 'vcf' block"
        )

    try:
        chrom = str(vcf["chr"])
        pos = int(vcf["pos"])
        ref = str(vcf["ref"])
        alt = str(vcf["alt"])
    except (KeyError, ValueError, TypeError) as exc:
        raise HGVSResolutionError(
            f"VariantValidator 'vcf' block is missing or malformed fields: {vcf!r}"
        ) from exc

    return chrom, pos, ref, alt


def resolve_coordinates(
    hgvs_c_with_transcript: str,
    genome_build: str = "GRCh38",
    base_url: str = DEFAULT_BASE_URL,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> ResolvedCoordinates:
    """Resolve a transcript-qualified HGVS coding variant description to
    strand-corrected genomic coordinates via the VariantValidator REST API.

    Args:
        hgvs_c_with_transcript: e.g. "NM_000162.5:c.1021G>A". Must include
            the transcript accession -- this function does not select one
            for you.
        genome_build: e.g. "GRCh38" (default) or "GRCh37".
        base_url: VariantValidator REST API base URL. Overridable for
            testing or for pointing at a self-hosted instance.
        timeout_seconds: HTTP request timeout.
        session: optional requests.Session to reuse connections /
            simplify mocking in tests. A new session is created if not
            supplied.
        debug: if True, print the raw response's top-level keys and any
            validation_warnings when something goes wrong, instead of
            just the summary error message.

    Returns:
        ResolvedCoordinates with strand-corrected chrom/pos/ref/alt.

    Raises:
        HGVSResolutionError: on any failure to obtain and verify a
            confident result -- malformed input, network/HTTP failure,
            unexpected response shape, or a response that fails the
            stale-cache verification check. Callers should catch this
            and fall back to manually-provided coordinates.
    """
    transcript, hgvs_c = _parse_transcript_and_hgvs(hgvs_c_with_transcript)
    genome_build = _normalize_genome_build(genome_build)

    http = session or requests.Session()

    url = (
        f"{base_url.rstrip('/')}/VariantValidator/variantvalidator/"
        f"{quote(genome_build, safe='')}/"
        f"{quote(hgvs_c_with_transcript, safe='')}/mane"
    )

    if debug:
        print(f"[hgvs_source] GET {url}")

    try:
        response = http.get(url, timeout=timeout_seconds)
    except requests.RequestException as exc:
        raise HGVSResolutionError(
            f"Network error calling VariantValidator for {hgvs_c_with_transcript!r}: {exc}"
        ) from exc

    if debug:
        print(f"[hgvs_source] status {response.status_code}")

    if response.status_code != 200:
        raise HGVSResolutionError(
            f"VariantValidator returned HTTP {response.status_code} for "
            f"{hgvs_c_with_transcript!r}: {response.text[:500]}"
        )

    try:
        payload = response.json()
    except ValueError as exc:
        raise HGVSResolutionError(
            f"VariantValidator response was not valid JSON for {hgvs_c_with_transcript!r}"
        ) from exc

    if debug:
        print(f"[hgvs_source] response top-level keys: {list(payload.keys())!r}")

    entry = _find_matching_result(payload, transcript)

    if debug:
        print(f"[hgvs_source] matched entry keys: {list(entry.keys())!r}")
        print(f"[hgvs_source] entry validation_warnings: {entry.get('validation_warnings')!r}")

    validated_description = _verify_match(entry, transcript, hgvs_c)
    chrom, pos, ref, alt = _extract_vcf_fields(entry, genome_build, debug=debug)

    return ResolvedCoordinates(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        genome_build=genome_build,
        transcript=transcript,
        hgvs_c=hgvs_c,
        validated_description=validated_description,
    )


# ---------------------------------------------------------------------------
# Gene symbol -> MANE Select transcript resolution
# ---------------------------------------------------------------------------
#
# NOTE ON LIVE VERIFICATION: this sandbox has no outbound network access, so
# the parsing below was written against VariantValidator's documented
# gene2transcripts response shape, not a live call. VariantValidator's
# response nests transcripts under a "transcripts" list, with a MANE Select
# flag living somewhere under each transcript's "annotations" block (exact
# key/casing has shifted across API versions in the past, similar to the
# dbNSFP field-path drift noted in myvariant_source.py). _find_mane_select()
# below therefore checks several plausible flag locations rather than one,
# but -- same caveat as clinvar_source.py's XPath comments -- confirm this
# against one live response for your gene of interest before trusting it on
# a real case, and adjust the checked keys if none of them match.

GENE2TRANSCRIPTS_PATH = "VariantValidator/tools/gene2transcripts"

_MANE_SELECT_LABELS = {"mane select", "mane_select"}


def _entry_is_mane_select(entry: dict[str, Any]) -> bool:
    """Check the plausible locations for a MANE Select flag on a single
    transcript entry from a gene2transcripts response.
    """
    annotations = entry.get("annotations") if isinstance(entry.get("annotations"), dict) else {}

    # Boolean-flag style: {"annotations": {"mane_select": true}}
    if annotations.get("mane_select") is True:
        return True
    if entry.get("mane_select") is True:
        return True

    # Status-string style: {"annotations": {"mane_status": "MANE Select"}}
    # or {"mane_status": "MANE Select"} / {"mane": "select"}
    for key, container in (
        ("mane_status", annotations),
        ("mane", annotations),
        ("mane_status", entry),
        ("mane", entry),
    ):
        value = container.get(key)
        if isinstance(value, str) and value.strip().lower() in _MANE_SELECT_LABELS:
            return True

    return False


def _find_mane_select_transcript(payload: dict[str, Any], gene_symbol: str) -> str:
    """Extract the MANE Select transcript accession from a
    gene2transcripts response payload.

    Raises HGVSResolutionError if the payload doesn't contain exactly one
    transcript flagged as MANE Select -- we do not guess among ambiguous
    candidates.
    """
    if not isinstance(payload, dict):
        raise HGVSResolutionError(
            f"Unexpected gene2transcripts response shape for {gene_symbol!r} "
            f"(not a dict): {type(payload)}"
        )

    # Response may be a single gene object, or a list of gene objects
    # (VariantValidator can return a list even for one query symbol).
    gene_entries: list[dict[str, Any]]
    if "transcripts" in payload:
        gene_entries = [payload]
    elif isinstance(payload.get("data"), list):
        gene_entries = [g for g in payload["data"] if isinstance(g, dict)]
    else:
        # Last resort: maybe the payload itself is shaped as {symbol: {...}}
        gene_entries = [v for v in payload.values() if isinstance(v, dict) and "transcripts" in v]

    if not gene_entries:
        raise HGVSResolutionError(
            f"Could not find a 'transcripts' list in gene2transcripts response "
            f"for gene {gene_symbol!r}. Response keys: {list(payload.keys())!r}"
        )

    mane_accessions: list[str] = []
    for gene_entry in gene_entries:
        transcripts = gene_entry.get("transcripts")
        if not isinstance(transcripts, list):
            continue
        for entry in transcripts:
            if not isinstance(entry, dict):
                continue
            if _entry_is_mane_select(entry):
                accession = entry.get("reference") or entry.get("transcript") or entry.get("accession")
                if accession:
                    mane_accessions.append(str(accession))

    unique_accessions = sorted(set(mane_accessions))

    if len(unique_accessions) == 1:
        return unique_accessions[0]

    if not unique_accessions:
        raise HGVSResolutionError(
            f"No transcript flagged as MANE Select was found for gene "
            f"{gene_symbol!r}. This can mean the gene has no MANE Select "
            f"transcript yet, the gene symbol wasn't recognised, or the "
            f"flag lives at a response path this module doesn't check yet "
            f"(see the NOTE ON LIVE VERIFICATION above _find_mane_select_transcript). "
            f"Supply the transcript accession explicitly instead, e.g. via "
            f"resolve_coordinates(\"NM_XXXXXX.X:c...\", ...)."
        )

    raise HGVSResolutionError(
        f"Found multiple different transcripts flagged as MANE Select for "
        f"gene {gene_symbol!r}: {unique_accessions!r}. Refusing to guess -- "
        f"supply the transcript accession explicitly instead."
    )


def resolve_mane_transcript(
    gene_symbol: str,
    base_url: str = DEFAULT_BASE_URL,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> str:
    """Look up the MANE Select RefSeq transcript accession for a gene
    symbol via VariantValidator's gene2transcripts endpoint.

    Args:
        gene_symbol: e.g. "GCK", "BRCA2".
        base_url, timeout_seconds, session: as in resolve_coordinates().
        debug: if True, print the resolved transcript accession, and (on
            failure) the raw response's top-level shape.

    Returns:
        The transcript accession string, e.g. "NM_000162.5".

    Raises:
        HGVSResolutionError: on network/HTTP failure, unexpected response
            shape, or if a single MANE Select transcript can't be
            confidently identified.
    """
    http = session or requests.Session()
    url = f"{base_url.rstrip('/')}/{GENE2TRANSCRIPTS_PATH}/{quote(gene_symbol, safe='')}"

    if debug:
        print(f"[hgvs_source] GET {url}")

    try:
        response = http.get(url, timeout=timeout_seconds)
    except requests.RequestException as exc:
        raise HGVSResolutionError(
            f"Network error calling VariantValidator gene2transcripts for "
            f"{gene_symbol!r}: {exc}"
        ) from exc

    if response.status_code != 200:
        raise HGVSResolutionError(
            f"VariantValidator gene2transcripts returned HTTP "
            f"{response.status_code} for {gene_symbol!r}: {response.text[:500]}"
        )

    try:
        payload = response.json()
    except ValueError as exc:
        raise HGVSResolutionError(
            f"VariantValidator gene2transcripts response was not valid JSON "
            f"for {gene_symbol!r}"
        ) from exc

    if debug:
        print(f"[hgvs_source] gene2transcripts response top-level keys: {list(payload.keys())!r}")

    transcript = _find_mane_select_transcript(payload, gene_symbol)

    if debug:
        print(f"[hgvs_source] MANE Select transcript for {gene_symbol!r}: {transcript}")

    return transcript


def resolve_coordinates_from_gene(
    gene_symbol: str,
    hgvs_c: str,
    genome_build: str = "GRCh38",
    base_url: str = DEFAULT_BASE_URL,
    timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    session: Optional["requests.Session"] = None,
    debug: bool = False,
) -> ResolvedCoordinates:
    """Resolve genomic coordinates from just a gene symbol and c. notation
    (e.g. gene_symbol="GCK", hgvs_c="c.1021G>A"), with no transcript
    accession required from the caller.

    This is a convenience wrapper around resolve_mane_transcript() +
    resolve_coordinates(): it looks up the gene's MANE Select transcript,
    then resolves coordinates against it. The transcript actually used is
    always available on the result as `.transcript`, so it can be shown
    to a reviewer rather than staying hidden.

    Prefer resolve_coordinates() with an explicit transcript accession
    when you already have one -- it's one fewer network round trip and
    removes any ambiguity about transcript selection.

    Raises:
        HGVSResolutionError: if the gene's MANE Select transcript can't
            be resolved, or if coordinate resolution against it fails
            (including the stale-cache verification in
            resolve_coordinates()). There is no coordinate fallback here
            -- gene symbol + c. notation is the only input -- so callers
            should let this propagate rather than swallow it.
    """
    transcript = resolve_mane_transcript(
        gene_symbol, base_url=base_url, timeout_seconds=timeout_seconds, session=session,
        debug=debug,
    )
    hgvs_c_with_transcript = f"{transcript}:{hgvs_c.strip()}"
    return resolve_coordinates(
        hgvs_c_with_transcript,
        genome_build=genome_build,
        base_url=base_url,
        timeout_seconds=timeout_seconds,
        session=session,
        debug=debug,
    )
