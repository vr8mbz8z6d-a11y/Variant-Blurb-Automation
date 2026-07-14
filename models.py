"""
Data schema for the variant blurb pipeline.

Every lookup source (gnomAD, ClinVar, PubMed, REVEL, SpliceAI) fills in a
piece of this record. Fields default to None so that a sentence module can
cleanly skip a fact that isn't available, rather than guessing or rendering
a blank.
"""
from __future__ import annotations
from typing import List, Optional
from pydantic import BaseModel, Field


class ClinVarSubmission(BaseModel):
    """A single lab/submitter's reported classification for this variant."""
    submitter: str
    classification: str
    review_status: Optional[str] = None


class GnomadData(BaseModel):
    af: Optional[float] = None          # allele frequency, e.g. 0.0000124
    allele_count: Optional[int] = None  # number of observed alt alleles
    allele_number: Optional[int] = None # total alleles genotyped at this site
    homozygote_count: Optional[int] = None


class ClinVarData(BaseModel):
    variation_id: Optional[str] = None              # ClinVar Variation/VCV ID
    aggregate_classification: Optional[str] = None    # ClinVar's overall call
    review_status: Optional[str] = None
    condition: Optional[str] = None                   # associated disease/condition name
    status: Optional[str] = None                       # e.g. "not_found" -- set by clinvar_source.py
    submissions: List[ClinVarSubmission] = Field(default_factory=list)


class RevelData(BaseModel):
    score: Optional[float] = None  # 0-1, only meaningful for missense variants


class SpliceAIData(BaseModel):
    """Only populated when variant_type == 'splice' (or splice-relevant)."""
    donor_loss: Optional[float] = None
    donor_gain: Optional[float] = None
    acceptor_loss: Optional[float] = None
    acceptor_gain: Optional[float] = None

    @property
    def max_delta(self) -> Optional[float]:
        vals = [v for v in (self.donor_loss, self.donor_gain,
                             self.acceptor_loss, self.acceptor_gain) if v is not None]
        return max(vals) if vals else None


class SpliceSiteContext(BaseModel):
    """
    Describes WHICH exon boundary a splice-region variant sits at, for a
    sentence like "This variant alters the invariant +2 nucleotide of the
    splice donor site immediately following exon 12." Only populated for
    variant_type == 'splice'. Derived from Ensembl VEP's exon/intron
    numbering (see ensembl_hgvs_source.get_splice_site_context) -- left
    None if that data isn't available, rather than guessing an exon
    number.
    """
    offset: Optional[str] = None       # e.g. "+2" or "-1", straight from the HGVS c. notation
    site_type: Optional[str] = None    # "donor" or "acceptor"
    nearest_exon: Optional[int] = None  # the exon number adjacent to the splice site


class NMDPrediction(BaseModel):
    """
    Predicted nonsense-mediated decay (NMD) outcome for a premature
    termination codon (nonsense/frameshift variants), per the standard
    Nagy & Maquat rule used in clinical variant interpretation (and the
    basis for ACMG PVS1 guidance):
      - PTC in the LAST exon -> escapes NMD
      - PTC in the second-to-last exon, within the last ~50 bp of it -> escapes NMD
      - Everywhere else -> triggers NMD
    Only populated for variant_type in ('nonsense', 'frameshift'). Left
    with escapes_nmd=None if the underlying exon-structure/position data
    couldn't be resolved -- never guesses.
    """
    escapes_nmd: Optional[bool] = None  # True = predicted to escape NMD, False = predicted to trigger it
    reason: Optional[str] = None        # e.g. "last exon", "within 50 bp of the final exon-exon junction"
    ptc_exon_number: Optional[int] = None
    total_exon_count: Optional[int] = None
    distance_to_last_junction: Optional[int] = None  # nt; None if PTC is in the last exon (not applicable)


class CoLocatedVariant(BaseModel):
    """
    A DIFFERENT variant affecting the same amino acid position as the
    variant being reported (e.g. query is p.Asn1355Lysfs*10; this
    represents a separately-reported p.Asn1355Lys at the same residue).
    Confirmed against your own reference blurb: "Another variant in the
    same position, p.Arg5105Gln, has conflicting classifications of
    pathogenicity by other clinical laboratories (ClinVar ID: 69441)."
    """
    protein_change: str          # e.g. "p.Asn1355Lys"
    hgvs_c: Optional[str] = None  # e.g. "c.8379A>C" -- CONFIRMED BUG FIX: needed
                                   # to disambiguate co-located variants whose
                                   # protein notation collapses to an identical
                                   # string (most notably synonymous variants,
                                   # where c.8379A>C, c.8379A>T, and c.8379A>G
                                   # all display as the SAME "p.Gly2793=" even
                                   # though they're three distinct variants --
                                   # confirmed against a real live case, BRCA2
                                   # p.Gly2793Arg's co-located search).
    classification: Optional[str] = None
    clinvar_id: Optional[str] = None


class FunctionalStudyFlag(BaseModel):
    """
    Deliberately minimal placeholder. Functional-study extraction is being
    built separately (see prior discussion) and is NOT auto-populated by
    this pipeline. Left here so the schema/template layer doesn't need to
    change shape later when that piece is added.
    """
    pmids: List[str] = Field(default_factory=list)
    note: str = "Functional study detection not yet wired into this pipeline."


class VariantRecord(BaseModel):
    """The single object that flows through lookup -> template rendering."""
    # input identifiers
    chrom: str
    pos: int
    ref: str
    alt: str
    genome_build: str = "hg38"
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    hgvs_transcript: Optional[str] = None  # transcript accession actually resolved against, e.g. "NM_000162.5"
    gene_strand: Optional[int] = None  # +1 or -1 from Ensembl VEP, used for strand-aware splice-site wording
    gene: Optional[str] = None
    variant_type: Optional[str] = None  # missense | nonsense | frameshift | splice | splice site variant | synonymous | indel | unknown
    clinical_classification: Optional[str] = None  # user-entered final classification, e.g. Pathogenic, Likely benign, VUS

    # populated by source modules
    gnomad: GnomadData = Field(default_factory=GnomadData)
    clinvar: ClinVarData = Field(default_factory=ClinVarData)
    revel: RevelData = Field(default_factory=RevelData)
    spliceai: SpliceAIData = Field(default_factory=SpliceAIData)
    splice_context: SpliceSiteContext = Field(default_factory=SpliceSiteContext)
    nmd: NMDPrediction = Field(default_factory=NMDPrediction)
    co_located_variants: List[CoLocatedVariant] = Field(default_factory=list)
    pubmed_ids: List[str] = Field(default_factory=list)
    functional_study: FunctionalStudyFlag = Field(default_factory=FunctionalStudyFlag)

    def myvariant_hgvs_id(self) -> str:
        """chr1:g.12345A>T style ID, the primary key MyVariant.info expects."""
        chrom = self.chrom if self.chrom.lower().startswith("chr") else f"chr{self.chrom}"
        return f"{chrom}:g.{self.pos}{self.ref}>{self.alt}"
