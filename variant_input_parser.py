"""
Parses a single, ClinVar-style combined variant string -- the format
ClinVar itself displays, e.g.:
    "NM_000326.5(RLBP1):c.753C>A; p.(Tyr251*)"
into its component pieces (transcript, gene, hgvs_c, hgvs_p), so you can
paste one string instead of re-typing gene/transcript/hgvs_c/hgvs_p
separately.

Tested against the real format variations actually seen in this project:
  - "NM_000326.5(RLBP1):c.753C>A; p.(Tyr251*)"      (semicolon + parenthesized p.)
  - "NM_022356.4(P3H1):c.1838+2T>C"                  (splice variant, no p. at all)
  - "NM_000059.4(BRCA2):c.8377G>A (p.Gly2793Arg)"    (space + unparenthesized p.)
  - "NM_022356.4:c.1838+2T>C"                        (no gene symbol given)
  - "NM_000326.5(RLBP1):c.753C>A p.Tyr251*"          (no punctuation before p.)

A transcript accession is always REQUIRED -- this parser never guesses
one, consistent with the rest of this pipeline. Raises ValueError with a
clear message (showing the expected format) on anything it can't parse,
rather than silently returning partial/wrong data.

The protein notation's own parentheses (the HGVS convention for a
"predicted, not directly confirmed" consequence, e.g. "p.(Tyr251*)") are
stripped for internal storage -- every reference blurb example built
throughout this project uses the non-parenthetical style in the
rendered sentence ("The p.Gly2793Arg missense variant..."), and keeping
the parens would also break existing regex-based variant-type
classification in variant_type.py, which expects "p." to be immediately
followed by the amino-acid letters.
"""
from __future__ import annotations
import re
from dataclasses import dataclass
from typing import Optional


class VariantStringParseError(ValueError):
    """Raised when a combined variant string can't be confidently parsed."""


@dataclass
class ParsedVariantInput:
    transcript: str
    gene: Optional[str]
    hgvs_c: str
    hgvs_p: Optional[str]

    @property
    def hgvs_c_with_transcript(self) -> str:
        return f"{self.transcript}:{self.hgvs_c}"


_COMBINED_PATTERN = re.compile(
    r"^\s*(?P<transcript>[A-Z]{2}_\d+\.\d+)"           # NM_000326.5
    r"(?:\((?P<gene>[A-Za-z0-9\-]+)\))?"                 # (RLBP1) -- optional
    r"\s*:\s*"                                            # :
    r"(?P<hgvs_c>c\.[^\s;(:]+)"                           # c.753C>A -- FIXED: also stop at ':'
    r"(?:[\s;,:]*\(?p\.\(?(?P<hgvs_p_inner>[^)\s]+)\)?\)?)?"  # ; p.(Tyr251*) or : p.Tyr251* -- optional
    r"\s*$"
)


def parse_combined_variant_string(raw: str) -> ParsedVariantInput:
    """
    Parse a combined ClinVar-style variant string into its parts.

    Raises:
        VariantStringParseError: if the string doesn't contain at least
            a transcript accession and a c. notation change -- these two
            are always required. Gene symbol and protein notation are
            optional (e.g. splice variants typically have no p. change).
    """
    if not raw or not raw.strip():
        raise VariantStringParseError("No variant string provided.")

    match = _COMBINED_PATTERN.match(raw.strip())
    if not match:
        raise VariantStringParseError(
            f"Could not parse {raw!r} as a variant description. Expected a "
            f"format like 'NM_000326.5(RLBP1):c.753C>A; p.(Tyr251*)' -- a "
            f"transcript accession (e.g. NM_000326.5), optionally followed by "
            f"the gene symbol in parentheses, a colon, the c. notation change, "
            f"and optionally a p. notation change. A transcript accession and "
            f"c. notation are always required; this tool never guesses either."
        )

    hgvs_p = f"p.{match.group('hgvs_p_inner')}" if match.group("hgvs_p_inner") else None

    return ParsedVariantInput(
        transcript=match.group("transcript"),
        gene=match.group("gene"),
        hgvs_c=match.group("hgvs_c"),
        hgvs_p=hgvs_p,
    )
