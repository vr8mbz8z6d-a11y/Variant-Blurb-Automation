from __future__ import annotations
import re
from models import VariantRecord


def nonsense_mechanism_sentence(v: VariantRecord) -> str | None:
    if v.variant_type != "nonsense":
        return None

    if not v.hgvs_p:
        return (
            "This nonsense variant leads to a premature termination codon, "
            "which is predicted to lead to a truncated or absent protein."
        )

    # Extract the amino acid position from p.Gln3982*, p.Arg97Ter, etc.
    match = re.search(r"[A-Za-z]+(\d+)(?:\*|Ter)", v.hgvs_p)

    if match:
        position = match.group(1)
        return (
            f"This nonsense variant leads to a premature termination codon "
            f"at position {position}, which is predicted to lead to a "
            f"truncated or absent protein."
        )

    # Fallback if parsing fails
    return (
        "This nonsense variant leads to a premature termination codon, "
        "which is predicted to lead to a truncated or absent protein."
    )

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


def _clinvar_own_label(v: VariantRecord) -> str | None:
    """
    Determine the query variant's OWN ClinVar-derived classification
    label, using the aggregate call if present, else a compact label
    built from distinct per-submission classifications. Extracted as a
    shared helper so clinvar_other_labs_sentence and
    co_located_variants_sentence (which folds this same label into a
    combined statement when co-located variants exist) can't disagree
    with each other about what the query variant's own label is.
    """
    label = (v.clinvar.aggregate_classification or "").strip()
    if label:
        return label
    seen = []
    for s in v.clinvar.submissions:
        c = (s.classification or "").strip().capitalize()
        if c and c not in seen:  # skip blank/whitespace-only entries
            seen.append(c)
    return "/ ".join(seen) if seen else None


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
    
    label = _clinvar_own_label(v)

    # If co-located variants were found, the combined
    # co_located_variants_sentence takes over stating this variant's own
    # classification (folded in alongside the co-located ones, e.g.
    # "This and other variants at the same codon... have been classified
    # as X/Y..."), so this sentence steps aside to avoid stating the same
    # classification twice in one blurb.
    if label and v.co_located_variants:
        return None

    if not label:
        return fallback_sentence

    return f"This variant has been classified as {label} by other clinical laboratories{id_str}."


def _join_natural(items: list[str]) -> str:
    """Join strings in natural English list style: 'A' / 'A and B' /
    'A, B, and C'. Used for both the protein-change list and the
    ClinVar ID list in co_located_variants_sentence."""
    items = list(items)
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return f"{items[0]} and {items[1]}"
    return ", ".join(items[:-1]) + f", and {items[-1]}"


def co_located_variants_sentence(v: VariantRecord) -> str | None:
    """
    Reports OTHER ClinVar variants affecting the SAME codon as this one,
    COMBINED with this variant's own classification into one statement,
    e.g.:
      "This and other variants at the same codon (p.Arg93Cys and
       Arg93His) have been classified as Pathogenic/Likely pathogenic
       by other clinical laboratories (ClinVar ID: 553695 and 1880)."

    Confirmed against a real, explicit target example (MMUT p.Arg93Ser,
    query classification Pathogenic; co-located Arg93Cys=Likely
    pathogenic, Arg93His=Pathogenic) -- deduped, query-first ordering of
    distinct classifications reproduces "Pathogenic/Likely pathogenic"
    exactly (query's own "Pathogenic" first, then Cys's "Likely
    pathogenic" as the next new distinct value; His's "Pathogenic" is
    already seen, so it isn't repeated).

    IMPORTANT: this REPLACES clinvar_other_labs_sentence's Case 3 output
    when co-located variants are found (see that function's own gate) --
    it folds the query variant's own classification in here instead of
    stating it twice. The ClinVar ID list includes the query variant's
    own ID first, followed by each co-located variant's -- e.g. for the
    MMUT example, "(ClinVar ID: 829883, 553695, and 1880)" where 829883
    is the query's own record. The protein-change list, by contrast,
    does NOT repeat the query's own protein change (only the co-located
    variants' changes appear there), since "This" already refers to it
    and restating "p.Arg93Ser and Arg93Cys and Arg93His" would be
    redundant with the blurb's own subject line.

    Grammar: "This and another variant" (singular) for exactly one
    co-located variant found; "This and other variants" (plural) for
    more than one. The verb is always plural ("have"), since "this AND
    [at least one other]" is always a 2+-item subject either way.

    Falls back to "have also been reported in ClinVar" (no
    classification claimed) if NEITHER the query nor any co-located
    variant has a known classification -- never asserts a label that
    isn't actually known for anyone in the group.

    Omitted entirely if no co-located variants were found (the common
    case) -- this data is populated by find_co_located_variants in
    clinvar_source.py.
    """
    if not v.co_located_variants:
        return None

    own_label = _clinvar_own_label(v)

    # Distinct classifications across the whole group, query's own
    # first, then each co-located variant's in the order found --
    # deduped so an identical label isn't repeated.
    combined_labels: list[str] = []
    if own_label:
        combined_labels.append(own_label)
    for cv in v.co_located_variants:
        if cv.classification and cv.classification not in combined_labels:
            combined_labels.append(cv.classification)

    # Protein-change list: keep "p." on the first item only, matching
    # natural English list convention (e.g. "p.Arg93Cys and Arg93His",
    # not "p.Arg93Cys and p.Arg93His").
    #
    # CONFIRMED BUG FIX: protein notation alone can be genuinely
    # ambiguous -- most notably for synonymous variants, where
    # c.8379A>C, c.8379A>T, and c.8379A>G all display as the SAME
    # "p.Gly2793=" despite being three distinct variants (confirmed
    # against a real live case, BRCA2 p.Gly2793Arg's co-located search,
    # which rendered three indistinguishable "Gly2793=" entries). Any
    # co-located variant whose protein_change string is shared by
    # another co-located variant IN THIS SAME GROUP falls back to
    # showing its c. notation instead (e.g. "c.8379A>C" rather than
    # "Gly2793="), so each entry stays individually identifiable. A
    # protein_change that's unique within the group is displayed as
    # before -- this only changes behavior for the genuinely ambiguous
    # case, not the common case your original examples were built on.
    protein_change_counts: dict[str, int] = {}
    for cv in v.co_located_variants:
        protein_change_counts[cv.protein_change] = protein_change_counts.get(cv.protein_change, 0) + 1

    raw_changes = []
    for cv in v.co_located_variants:
        if protein_change_counts[cv.protein_change] > 1 and cv.hgvs_c:
            raw_changes.append(cv.hgvs_c)  # disambiguate with c. notation
        else:
            raw_changes.append(cv.protein_change)

    display_changes = [raw_changes[0]] + [
        pc[2:] if pc.startswith("p.") else pc for pc in raw_changes[1:]
    ]
    protein_list_str = _join_natural(display_changes)

    # CONFIRMED FIX: "This and other variants..." falsely implies the
    # query variant itself is part of the classified group. When the
    # query has no own classification to contribute (own_label is None
    # -- either genuinely absent from ClinVar, or present but with no
    # submitted interpretations), drop "This and" entirely and start
    # with just "Other variants"/"Another variant" instead, since only
    # the CO-LOCATED variants are actually the ones being described as
    # classified. This also resolves a real contradiction: without this,
    # the blurb could say "This variant has not been reported in
    # ClinVar." immediately followed by "This and other variants...
    # have been classified as X" -- directly conflicting statements
    # about the same "This".
    if own_label:
        subject = "This and another variant" if len(v.co_located_variants) == 1 else "This and other variants"
    else:
        subject = "Another variant" if len(v.co_located_variants) == 1 else "Other variants"

    # Subject-verb agreement: only the bare "Another variant" (singular,
    # no "This and") takes "has" -- every other case ("This and another
    # variant" = 2 items, "This and other variants" / "Other variants" =
    # plural) takes "have".
    verb = "has" if subject == "Another variant" else "have"

    # Query's own ClinVar ID first (matching "This" being mentioned
    # first in the subject), then each co-located variant's ID in order.
    # Gated on own_label (not just presence of a ClinVar record): the
    # query's own ID is only included when it's actually part of the
    # "This and..." group being described -- consistent with the
    # subject-phrasing fix above. If the query is in ClinVar but has no
    # classification (Case 2-style), its ID would otherwise appear
    # alongside a sentence that no longer even mentions "This".
    all_ids = ([v.clinvar.variation_id] if own_label and v.clinvar.variation_id else []) + \
              [cv.clinvar_id for cv in v.co_located_variants if cv.clinvar_id]
    id_list_str = _join_natural(all_ids)
    id_clause = f" (ClinVar ID: {id_list_str})" if id_list_str else ""

    if combined_labels:
        classification_clause = (
            f"{verb} been classified as {'/'.join(combined_labels)} by other clinical laboratories"
        )
    else:
        classification_clause = f"{verb} also been reported in ClinVar"

    return f"{subject} at the same codon ({protein_list_str}) {classification_clause}{id_clause}."


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
    nonsense_mechanism_sentence,
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
