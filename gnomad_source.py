"""
gnomAD allele-frequency lookup via the official gnomAD GraphQL API (v4).

WHY THIS EXISTS (and why we don't just use MyVariant.info for gnomAD):
MyVariant.info's gnomad_exome / gnomad_genome fields come from gnomAD v2.1,
whose GENOME data was mapped to hg19 only. On the hg38 index you therefore
get the exome block but no matching genome block, so a "combined total"
can't be assembled there -- that's why the pipeline was reporting an
exome-only denominator (e.g. 2/251,416) instead of the true combined
figure.

UPDATED: this module used to compute the combined total itself by
manually summing the `exome` and `genome` blocks' ac/an/homozygote
counts. gnomAD's own GraphQL API actually exposes a dedicated `joint`
field that IS the authoritative combined figure -- confirmed against a
real, independently-verified working example (a third-party tool built
directly on gnomad-browser's GraphQL API) whose sample output's joint
AC/AN exactly equals the sum of its exome + genome AC/AN (395,737 +
52,299 = 448,036 joint AC; 1,315,438 + 151,162 = 1,466,600 joint AN).
Querying `joint` directly is preferred over manually summing ourselves:
it's the API's own authoritative calculation (matching gnomAD v4.1's
documented joint-AN methodology, which retains allele-number information
across all called sites, a subtlety a naive manual sum could miss for
some variants) and it's simpler and less error-prone than re-deriving it
client-side. Manual summing of exome+genome is kept ONLY as a fallback
for the rare case a `joint` block isn't present.

Also corrected the homozygote field name from `ac_hom` to
`homozygote_count`, confirmed against the same verified working example
-- the two may or may not be interchangeable aliases in gnomAD's schema,
but using the field name from confirmed-working code removes any doubt.

Query shape and field names verified against a real, working third-party
script (github.com/sfbizzari/gnomADv4-Batch-tool-pythonAPI), not just
documentation:
  variant(variantId: "2-233760233-C-CAT", dataset: gnomad_r4) {
    exome  { ac an homozygote_count }
    genome { ac an homozygote_count }
    joint  { ac an homozygote_count }
  }
Note: variantId is 'chrom-pos-ref-alt' with NO 'chr' prefix; all three
blocks may be null if the variant is absent from that dataset/grouping.

RATE LIMITING: the public API blocks aggressive use (reports of blocks
after ~10 rapid queries). This module is for single-variant/interactive
use, with a delay between calls. For batch annotation, gnomAD recommends
downloading their Hail tables / VCFs instead.
"""
from __future__ import annotations
from typing import Optional
import time
import requests

from models import GnomadData

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
TIMEOUT_SECONDS = 20
MIN_SECONDS_BETWEEN_CALLS = 2  # be polite to the public endpoint

_VARIANT_QUERY = """
query VariantTotals($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    exome  { ac an homozygote_count }
    genome { ac an homozygote_count }
    joint  { ac an homozygote_count }
  }
}
"""


def _variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """gnomAD wants 'chrom-pos-ref-alt' with no 'chr' prefix."""
    chrom_clean = chrom.replace("chr", "").replace("Chr", "")
    return f"{chrom_clean}-{pos}-{ref}-{alt}"


def fetch_gnomad_variant(chrom: str, pos: int, ref: str, alt: str,
                          dataset: str = "gnomad_r4", debug: bool = False) -> tuple[Optional[dict], bool]:
    """
    Returns (variant_dict_or_None, confirmed_absent).

    confirmed_absent=True means gnomAD's API was successfully reached and
    explicitly reported this variant isn't in the dataset -- a real,
    reportable clinical fact ("absent from gnomAD"), not the same as a
    failed lookup. gnomAD signals this via a GraphQL "errors" array with
    message "Variant not found" alongside data.variant=null, even on a
    200 HTTP status -- confirmed against a live response:
    {"errors":[{"message":"Variant not found"}],"data":{"variant":null}}

    Earlier versions of this function collapsed "confirmed absent" and
    "lookup genuinely failed" into the same None return value, which
    meant a real, reportable absence silently produced NO sentence at
    all in the rendered blurb instead of the correct "not observed in
    gnomAD" statement.
    """
    variant_id = _variant_id(chrom, pos, ref, alt)
    payload = {"query": _VARIANT_QUERY,
               "variables": {"variantId": variant_id, "dataset": dataset}}
    try:
        resp = requests.post(GNOMAD_API_URL, json=payload, timeout=TIMEOUT_SECONDS)
        time.sleep(MIN_SECONDS_BETWEEN_CALLS)
        if debug:
            print(f"[gnomad_source] POST variant={variant_id} dataset={dataset} -> {resp.status_code}")
            print(f"[gnomad_source] response (first 600 chars): {resp.text[:600]}")
        if resp.status_code != 200:
            return None, False
        body = resp.json()
        errors = body.get("errors") or []
        if debug and errors:
            print(f"[gnomad_source] GraphQL errors: {errors}")
        variant = body.get("data", {}).get("variant")
        confirmed_absent = variant is None and any(
            "not found" in (e.get("message") or "").lower() for e in errors
        )
        return variant, confirmed_absent
    except (requests.RequestException, ValueError) as e:
        if debug:
            print(f"[gnomad_source] request failed: {e}")
        return None, False


def parse_combined_gnomad(variant: Optional[dict], confirmed_absent: bool = False,
                           debug: bool = False) -> GnomadData:
    """
    Return the combined (exome + genome) gnomAD figure, preferring the
    API's own authoritative `joint` block. Falls back to manually summing
    `exome` + `genome` only if `joint` isn't present in the response
    (e.g. an older dataset, or a schema variation not yet seen live).

    confirmed_absent=True (gnomAD explicitly reported "Variant not
    found", not a lookup failure) returns af=0.0 rather than an empty
    GnomadData -- this is what lets the "not observed in gnomAD" sentence
    render. Without this, a confirmed absence and a failed/unknown lookup
    were indistinguishable, and the real fact "checked, not present" was
    silently dropped from the blurb instead of being reported.
    """
    if confirmed_absent:
        return GnomadData(af=0.0, allele_count=0, allele_number=0, homozygote_count=0)

    if not variant:
        return GnomadData()

    joint = variant.get("joint")
    if joint and (joint.get("ac") is not None or joint.get("an") is not None):
        ac = joint.get("ac") or 0
        an = joint.get("an") or 0
        hom = joint.get("homozygote_count") or 0
        if debug:
            print(f"[gnomad_source] using authoritative 'joint' block: ac={ac}, an={an}, hom={hom}")
        af = (ac / an) if an > 0 else None
        return GnomadData(af=af, allele_count=ac, allele_number=an, homozygote_count=hom)

    # Fallback: no joint block present -- sum exome + genome ourselves.
    if debug:
        print("[gnomad_source] no 'joint' block in response -- falling back to "
              "manually summing exome + genome blocks")

    total_ac = 0
    total_an = 0
    total_hom = 0
    found_any = False

    for block_key in ("exome", "genome"):
        block = variant.get(block_key)
        if not block:
            continue
        ac = block.get("ac")
        an = block.get("an")
        hom = block.get("homozygote_count")
        if ac is None and an is None:
            continue
        found_any = True
        total_ac += ac or 0
        total_an += an or 0
        total_hom += hom or 0

    if not found_any:
        return GnomadData()

    af = (total_ac / total_an) if total_an > 0 else None
    return GnomadData(af=af, allele_count=total_ac, allele_number=total_an,
                      homozygote_count=total_hom)


def lookup_gnomad(chrom: str, pos: int, ref: str, alt: str,
                   genome_build: str = "hg38", debug: bool = False,
                   alternate_representations: tuple = ()) -> GnomadData:
    """
    Combined gnomAD v4 (GRCh38) allele frequency for a variant, using the
    API's own authoritative `joint` (exome+genome combined) figure.

    gnomAD v4 is GRCh38-only. If an hg19/GRCh37 build is passed we fall
    back to the v2.1.1 dataset (gnomad_r2_1), which is gnomAD's documented
    GRCh37 dataset -- but note v2 predates the `joint` field and this
    path is not the primary, tested one.

    alternate_representations: optional (chrom, pos, ref, alt) tuples
    describing the SAME variant at a different but equally valid VCF
    position -- see ResolvedCoordinates.alternate_representations. If the
    primary ID isn't found, each is tried before we conclude the variant
    is genuinely absent.

    WHY THE RETRY EXISTS: "absent from gnomAD" is a strong clinical claim
    that feeds directly into the rendered blurb, and for indels a
    not-found result is ambiguous between "genuinely never observed" and
    "observed, but indexed under a different left/right-shifted
    representation of the same change". CONFIRMED case: MSH6
    NM_000179.3:c.3261dup was reported as absent while gnomAD actually
    holds it (AC=91, AF=5.65e-05) under a left-aligned ID seven bases
    away. The primary fix for that is left-alignment upstream in
    ensembl_hgvs_source; this retry is the belt-and-braces second line,
    so a future normalization gap degrades into an extra API call rather
    than a false clinical statement.
    """
    dataset = "gnomad_r4"
    if genome_build.lower() in ("hg19", "grch37", "37"):
        dataset = "gnomad_r2_1"
        if debug:
            print("[gnomad_source] GRCh37 build requested -> using gnomad_r2_1 dataset "
                  "(predates the 'joint' field -- will fall back to manual summing)")
    variant, confirmed_absent = fetch_gnomad_variant(chrom, pos, ref, alt, dataset=dataset, debug=debug)

    if variant is None and alternate_representations:
        for alt_chrom, alt_pos, alt_ref, alt_alt in alternate_representations:
            if (alt_chrom, alt_pos, alt_ref, alt_alt) == (chrom, pos, ref, alt):
                continue
            if debug:
                print(f"[gnomad_source] primary ID {_variant_id(chrom, pos, ref, alt)} "
                      f"not found -- retrying equivalent representation "
                      f"{_variant_id(alt_chrom, alt_pos, alt_ref, alt_alt)} before "
                      f"concluding this variant is absent")
            retry_variant, retry_absent = fetch_gnomad_variant(
                alt_chrom, alt_pos, alt_ref, alt_alt, dataset=dataset, debug=debug)
            if retry_variant is not None:
                if debug:
                    print(f"[gnomad_source] FOUND under alternate representation "
                          f"{_variant_id(alt_chrom, alt_pos, alt_ref, alt_alt)} -- the "
                          f"primary coordinate was a normalization mismatch, NOT a real "
                          f"absence. Using the alternate's data.")
                variant, confirmed_absent = retry_variant, False
                break
            # Only keep asserting absence if every representation agrees.
            confirmed_absent = confirmed_absent and retry_absent

    if debug and confirmed_absent:
        print("[gnomad_source] gnomAD explicitly reports this variant as absent "
              "(not a lookup failure) -- will render as 'not observed in gnomAD'.")
    return parse_combined_gnomad(variant, confirmed_absent=confirmed_absent, debug=debug)
