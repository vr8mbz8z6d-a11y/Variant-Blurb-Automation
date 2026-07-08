"""
Lightweight variant-type classification from HGVS notation.

This is the regex-on-HGVS approach we discussed: fast, free, no genome
coordinates needed, and good enough for routing pipeline logic (e.g.
"only look up SpliceAI for splice variants"). It is NOT a replacement for
a real consequence predictor (VEP/SnpEff) -- see the caveats noted where
this is called.
"""
from __future__ import annotations
import re
from typing import Optional


def classify_variant_type(hgvs_p: Optional[str] = None, hgvs_c: Optional[str] = None) -> str:
    """
    Returns one of: missense, nonsense, frameshift, synonymous, indel,
    splice, splice site variant, unknown.

    Order of checks matters -- frameshift/nonsense markers are checked before
    falling back to a generic "two different residues = missense" assumption.
    """
    p = hgvs_p or ""
    c = hgvs_c or ""

    # Splice-region notation lives in the c. string, e.g. c.123+1G>A
    # or c.456-2A>G. Canonical +-1/+-2 variants stay in the existing
    # "splice" bucket; further intronic offsets are the separate
    # "splice site variant" type used for outside-consensus blurbs.
    splice_match = re.search(r"[+-](\d+)[ACGT]>[ACGT]", c)
    if splice_match:
        offset = int(splice_match.group(1))
        if offset > 2:
            return "splice site variant"
        return "splice"

    if "fs" in p:
        return "frameshift"

    if "*" in p or "Ter" in p:
        return "nonsense"

    if "del" in p or "ins" in p or "dup" in p:
        return "indel"

    if p.endswith("=") or "=" in p:
        return "synonymous"

    # p.Arg175His style: two different three-letter residues -> missense
    if re.search(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", p):
        ref_res = re.findall(r"p\.([A-Za-z]{3})", p)
        alt_res = re.findall(r"\d+([A-Za-z]{3})", p)
        if ref_res and alt_res and ref_res[0] != alt_res[0]:
            return "missense"

    return "unknown"
