"""
SpliceAI lookup, restricted to splice-relevant variants only (per your
criteria: "spliceAI should also be included but only for splice variants").

Uses the Broad Institute's public SpliceAI-lookup API
(https://github.com/broadinstitute/SpliceAI-lookup), the same engine behind
spliceailookup.broadinstitute.org. It's free and requires no API key, but
is explicitly rate-limited to "a handful of queries per user per minute" --
the README states it's intended for interactive/single-variant use, not
batch processing. For batch annotation of many variants, the maintainers'
own recommendation is to self-host the Docker image or run the SpliceAI
model directly; this module is written for the prototype/interactive case
only, with a deliberate delay between calls and a low retry count.

ENDPOINT CHANGE: previously pointed at the custom domain
spliceailookup-api.broadinstitute.org, which hit a hard TCP connect
timeout consistently across two independent networks (confirmed via curl
-v showing DNS resolved fine to 34.71.77.14 but the handshake itself
never completed), while the website itself (spliceailookup.broadinstitute.org)
worked fine in a browser. Per the project's own current README
(github.com/broadinstitute/SpliceAI-lookup), the actual, currently
documented API lives directly on Google Cloud Run's default domains,
with a SEPARATE hostname per genome build:
  https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=37&variant=...
  https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&variant=...
This is very likely what the live website's own JavaScript actually calls
-- which would explain the browser-vs-code discrepancy: the custom domain
may route through different infrastructure (a load balancer / domain
mapping) that was being filtered somewhere on the path, while the raw
Cloud Run domain wasn't.

NOT YET LIVE-VERIFIED: the response JSON shape from these specific
run.app URLs. The parsing below assumes the same shape as before
(a top-level "scores" list with DS_DL/DS_DG/DS_AL/DS_AG per transcript),
since this should be the same underlying Flask/Cloud Run app just served
from a different domain -- but this has not been confirmed against a
live response from this exact URL. If parsing comes back empty despite a
200 status, run with debug=True and check the raw payload structure
printed.
"""
from __future__ import annotations
from typing import Optional
import time
import requests

from models import SpliceAIData

# Separate base URL per genome build, per the current official README --
# NOT interchangeable via a query parameter alone; the hostname itself
# encodes the build.
SPLICEAI_BASE_URLS = {
    "37": "https://spliceai-37-xwkwwwxdwq-uc.a.run.app/spliceai/",
    "38": "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/",
}
TIMEOUT_SECONDS = 20
MIN_SECONDS_BETWEEN_CALLS = 3  # respect the documented per-minute rate limit


def lookup_spliceai(chrom: str, pos: int, ref: str, alt: str,
                     genome_build: str = "38", distance: int = 50,
                     debug: bool = False) -> SpliceAIData:
    """
    genome_build should be '37' or '38' -- selects BOTH the hostname and
    the `hg` query parameter, matching the exact documented usage (the
    README's own example includes `hg=38` in the query string even
    though the hostname already encodes build 38 -- kept as-is rather
    than "simplified", since redundant-looking parameters in someone
    else's API have turned out to matter before in this project).
    Variant format expected by the API: chrom-pos-ref-alt, e.g. chr8-140300616-T-G
    """
    base_url = SPLICEAI_BASE_URLS.get(str(genome_build))
    if base_url is None:
        if debug:
            print(f"[spliceai_source] unrecognized genome_build {genome_build!r}, "
                  f"expected '37' or '38' -- returning empty SpliceAIData")
        return SpliceAIData()

    chrom_clean = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    variant = f"{chrom_clean}-{pos}-{ref}-{alt}"
    params = {"hg": genome_build, "variant": variant, "distance": distance}

    if debug:
        print(f"[spliceai_source] GET {base_url} params={params}")

    try:
        resp = requests.get(base_url, params=params, timeout=TIMEOUT_SECONDS)
        time.sleep(MIN_SECONDS_BETWEEN_CALLS)
        if debug:
            print(f"[spliceai_source] status {resp.status_code}")
            print(f"[spliceai_source] response (first 500 chars): {resp.text[:500]}")
        if resp.status_code != 200:
            if debug:
                print(f"[spliceai_source] non-200 status -> returning empty SpliceAIData")
            return SpliceAIData()
        payload = resp.json()
    except (requests.RequestException, ValueError) as e:
        if debug:
            print(f"[spliceai_source] request failed: {e}")
        return SpliceAIData()

    # The API returns scores keyed per-transcript; we take the highest delta
    # per event type across transcripts, which is the conservative read.
    scores = payload.get("scores", [])
    if not scores:
        if debug:
            print(f"[spliceai_source] response had no 'scores' -- full payload: {payload!r}")
        return SpliceAIData()

    def _max_field(field: str) -> Optional[float]:
        # CONFIRMED BUG FIX: this API returns scores as STRINGS
        # (e.g. "0.08"), not JSON numbers, confirmed against a live
        # response: {"DS_DL": "0.08", ...}. The old isinstance(x, (int,
        # float)) check silently excluded every string value, so every
        # real score got discarded before max() ever ran, producing
        # None for all four fields despite valid data being present.
        vals = []
        for s in scores:
            raw = s.get(field)
            if raw is None:
                continue
            try:
                vals.append(float(raw))
            except (TypeError, ValueError):
                continue  # skip anything that isn't a parseable number
        return max(vals) if vals else None

    result = SpliceAIData(
        donor_loss=_max_field("DS_DL"),
        donor_gain=_max_field("DS_DG"),
        acceptor_loss=_max_field("DS_AL"),
        acceptor_gain=_max_field("DS_AG"),
    )
    if debug:
        print(f"[spliceai_source] parsed scores: donor_loss={result.donor_loss}, "
              f"donor_gain={result.donor_gain}, acceptor_loss={result.acceptor_loss}, "
              f"acceptor_gain={result.acceptor_gain}")
    return result
