"""
gnomAD allele-frequency and REVEL score lookup via MyVariant.info.

MyVariant.info (https://myvariant.info) is a free, public aggregation API
that merges ClinVar, gnomAD, dbNSFP (which includes REVEL), CADD, dbSNP and
~20 other sources by genomic HGVS ID. No API key is needed for normal use;
there's a documented rate limit (1,000 requests/IP/day for anonymous use).

We use it here for gnomAD + REVEL specifically because both already live in
one indexed, queryable place -- no multi-GB file download required, which
matches the "skip the heavy infrastructure for v1" plan we discussed.

IMPORTANT -- ASSEMBLY DEFAULT: MyVariant.info's variant IDs default to
hg19 coordinates unless you explicitly pass `assembly=hg38` as a query
parameter (https://docs.myvariant.info/en/latest/doc/data.html: "The
default reference genome assembly is always human hg19 in MyVariant.info").
Passing an hg38 position without that parameter doesn't error -- it just
silently looks up the wrong genomic location on hg19's map and returns
empty results, which is exactly what happened the first time this ran.

NOTE ON FIELD PATHS: dbNSFP nests its sub-scores under `dbnsfp.<tool>`.
The exact REVEL sub-path has shifted slightly across dbNSFP/MyVariant.info
schema versions (older: dbnsfp.revel.score, some releases: dbnsfp.revel_score).
This module tries both and logs which one matched, so you can confirm against
a live response and pin the right path for your dbNSFP version.
"""
from __future__ import annotations
from typing import Optional
import requests

from models import GnomadData, RevelData

MYVARIANT_BASE = "https://myvariant.info/v1/variant"
TIMEOUT_SECONDS = 10


def _get(d: dict, *path, default=None):
    """Safely walk a nested dict by a sequence of keys."""
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def fetch_myvariant_record(hgvs_genomic_id: str, genome_build: str = "hg38",
                            debug: bool = False) -> Optional[dict]:
    """
    hgvs_genomic_id example: 'chr13:g.32370447G>A'
    genome_build: 'hg38' (default here) or 'hg19'/'hg37'. Passed straight
    through to MyVariant.info's `assembly` parameter -- without this, an
    hg38 position gets looked up against hg19 coordinates instead (see
    module docstring), which is the bug that produced empty results.
    """
    url = f"{MYVARIANT_BASE}/{hgvs_genomic_id}"
    params = {}
    if genome_build.lower() in ("hg38", "grch38"):
        params["assembly"] = "hg38"
    try:
        resp = requests.get(url, params=params, timeout=TIMEOUT_SECONDS)
        if debug:
            print(f"[myvariant_source] GET {resp.url} -> status {resp.status_code}")
            print(f"[myvariant_source] response (first 500 chars): {resp.text[:500]}")
        if resp.status_code != 200:
            return None
        data = resp.json()
        # MyVariant.info returns a 404-style body for unknown variants
        if isinstance(data, dict) and data.get("success") is False:
            if debug:
                print("[myvariant_source] variant not found at this coordinate/assembly")
            return None
        return data
    except requests.RequestException as e:
        if debug:
            print(f"[myvariant_source] request failed: {e}")
        return None


def parse_gnomad(record: dict) -> GnomadData:
    """
    Combine gnomAD exome + genome data into a single total, matching what
    the gnomAD browser shows as the overall figure.

    MyVariant.info exposes 'gnomad_exome' and 'gnomad_genome' as two
    SEPARATE blocks, each with its own ac/an/hom. gnomAD's own "total"
    (and the denominator you'd cite in a report, e.g. 2/1,613,808) is the
    SUM of both datasets -- confirmed by gnomAD v4.1 docs and Nirvana's:
    "When merging the genomes and exomes, the allele counts and allele
    numbers will be summed across both of the data sets."

    The earlier version of this function returned only the first block it
    found, which under-reported the denominator (e.g. 251,416 exome-only
    instead of the combined total). We now sum whichever blocks are
    present. AF is recomputed as combined_ac / combined_an so it stays
    consistent with the summed counts, rather than using either dataset's
    standalone af.
    """
    total_ac = 0
    total_an = 0
    total_hom = 0
    found_any = False

    for key in ("gnomad_exome", "gnomad_genome"):
        block = record.get(key)
        if not block:
            continue
        ac = _get(block, "ac", "ac")
        an = _get(block, "an", "an")
        hom = _get(block, "hom", "hom")
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


def parse_revel(record: dict) -> RevelData:
    dbnsfp = record.get("dbnsfp", {})
    if not dbnsfp:
        return RevelData()

    score = _get(dbnsfp, "revel", "score")
    if score is None:
        score = dbnsfp.get("revel_score")
    if score is None:
        return RevelData()

    # dbNSFP sometimes returns a list when multiple transcripts overlap;
    # take the max as the more conservative (higher pathogenicity) estimate.
    if isinstance(score, list):
        score = max(score)
    return RevelData(score=float(score))


def lookup_gnomad_and_revel(hgvs_genomic_id: str, genome_build: str = "hg38",
                             debug: bool = False) -> tuple[GnomadData, RevelData]:
    record = fetch_myvariant_record(hgvs_genomic_id, genome_build=genome_build, debug=debug)
    if record is None:
        return GnomadData(), RevelData()
    return parse_gnomad(record), parse_revel(record)


def lookup_revel(hgvs_genomic_id: str, genome_build: str = "hg38",
                 debug: bool = False) -> RevelData:
    """
    REVEL-only lookup. gnomAD frequency now comes from the dedicated
    gnomAD GraphQL source (see gnomad_source.py) because MyVariant.info's
    hg38 index lacks the gnomAD genome block, so it can't give a true
    combined total. REVEL, however, lives in dbNSFP and is correctly
    available here, so we keep using MyVariant.info for it.
    """
    record = fetch_myvariant_record(hgvs_genomic_id, genome_build=genome_build, debug=debug)
    if record is None:
        return RevelData()
    return parse_revel(record)
