from pydantic import BaseModel, Field
from typing import Optional, Any
from tbp_parser.Utilities.gene_database import GeneDatabase
from tbp_parser.Utilities.helper import Helper
import logging

logger = logging.getLogger(__name__)

class Variant(BaseModel):
    """This class represents a single genetic drug-gene-variant combination that can be found in the Laboratorian report.

    Core attributes include fields expanded from the TBProfiler output JSON.

    This class assumes that the variant has already been expanded such that
    each Variant instance represents a single gene-drug association. These
    associations can be derived from the annotation and consequence field(s)
    in the TBProfiler output JSON.

    """
    # Core attributes from TBProfiler necessary to define a Variant
    sample_id: str
    pos: int
    depth: int
    freq: float
    gene_id: str
    gene_name: str
    type: str
    nucleotide_change: str
    protein_change: str
    drug: str
    confidence: str
    rationale: str = "NA"
    source: str = ""
    comment: str = ""

    # Derived fields (computed during init, excluded from serialization)
    read_support: Optional[float] = Field(default=None, exclude=True)
    gene_tier: Optional[str] = Field(default=None, exclude=True)
    absolute_start: Optional[int] = Field(default=None, exclude=True)
    absolute_end: Optional[int] = Field(default=None, exclude=True)

    # Fields set by VariantInterpreter (excluded from serialization)
    looker_interpretation: Optional[str] = Field(default=None, exclude=True)
    mdl_interpretation: Optional[str] = Field(default=None, exclude=True)

    # Fields set by VariantQC (excluded from serialization)
    fails_locus_qc: Optional[bool] = Field(default=None, exclude=True)
    fails_positional_qc: Optional[bool] = Field(default=None, exclude=True)
    warning: set[str] = Field(default_factory=set, exclude=True)

    model_config = {"extra": "ignore"}

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        # Calculate read_support if not provided
        if self.read_support is None:
            self.read_support = self.freq * self.depth

        # Derive computed attributes
        self.gene_tier = GeneDatabase.get_tier(self.gene_id)
        self.absolute_start, self.absolute_end = Helper.get_mutation_genomic_positions(self.pos, self.nucleotide_change)

        # Structural defaults
        if not self.protein_change:
            self.protein_change = "NA"
        if self.comment == "Not found in WHO catalogue":
            self.confidence = "No WHO annotation"

    def __str__(self) -> str:
        return f"Variant([{self.gene_name}][{self.gene_id}][{self.type}][{self.nucleotide_change}][{self.drug}][{self.confidence}])"

    def __repr__(self) -> str:
        return f"Variant([{self.gene_name}][{self.gene_id}][{self.type}][{self.nucleotide_change}][{self.drug}][{self.confidence}])"

    def __eq__(self, other):
        """Define equality based on some attributes"""
        if not isinstance(other, Variant):
            return False

        return (
            self.gene_id == other.gene_id and
            self.gene_name == other.gene_name and
            self.type == other.type and
            self.nucleotide_change == other.nucleotide_change and
            self.drug == other.drug
        )

    def __hash__(self):
        """Make the object hashable using some attributes"""
        return hash((
            self.gene_id,
            self.gene_name,
            self.type,
            self.nucleotide_change,
            self.drug
        ))

    @classmethod
    def from_thin_air(cls, sample_id: str, gene_id: str, gene_name: str, drug: str) -> 'Variant':
        """Create a Variant object from thin air.
        Will be converted to a WT/NA type Variant in VariantQC.

        Args:
            sample_id (str): The sample ID.
            gene_id (str): The locus tag of the gene.
            gene_name (str): The name of the gene.
            drug (str): The drug associated with the variant.

        Returns:
            Variant: A Variant object with appropriate fields set for WT/NA.
        """
        # model_construct will bypass validation and post-init processing,
        # allowing us to create a Variant with missing/NA fields that would normally fail validation
        return cls.model_construct(
            sample_id=sample_id,
            pos="NA",
            depth="NA",
            freq="NA",
            gene_id=gene_id,
            gene_name=gene_name,
            type="NA",
            nucleotide_change="NA",
            protein_change="NA",
            drug=drug,
            confidence="NA",
            read_support="NA",
        )

    def is_better_annotation_than(self, other: 'Variant') -> bool:
        annotation_rank = {
          "Assoc w R": 5,
          "Assoc w R - interim": 4,
          "Assoc w R - Interim": 4,
          "Uncertain significance": 3,
          "Not assoc w R - Interim": 2,
          "Not assoc w R - interim": 2,
          "Not assoc w R": 1,
          "Not found in WHO catalogue": 0,
          "No WHO annotation": -1, # given to synthetic variants
          "": -2, # empty confidence field
        }
        new_variant_rank = annotation_rank[self.confidence]
        existing_variant_rank = annotation_rank[other.confidence]

        # if it has been seen before, save the row with the more severe WHO confidence (higher value)
        if new_variant_rank > existing_variant_rank:
            logger.debug(f"{self} has a better annotation rank ({new_variant_rank}) than {other} ({existing_variant_rank})")
            return True

        # also preferentially keep WHO confidence rows over non-WHO confidence rows if the ranks are the same
        if new_variant_rank == existing_variant_rank:
            if "WHO" in self.source and "WHO" not in other.source:
                logger.debug(f"{self} has WHO source while {other} does not; preferring {self}")
                return True

        return False

    def has_better_read_support_than(self, other: 'Variant') -> bool:
        """Compare read support between two variants."""
        if self.read_support > other.read_support:  # type: ignore
            logger.debug(f"{self} has better read support ({self.read_support}) than {other} ({other.read_support})")
            return True

        return False

    def _is_synonymous(self) -> bool:
        """Check if mutation is synonymous."""
        return self.type == "synonymous_variant"

    def _is_upstream_gene_variant(self) -> bool:
        """Check if mutation is an upstream gene variant."""
        return "upstream_gene_variant" in self.type

    def _is_in_target_promoter(self, position_nt: list[int]) -> bool:
        """Check if variant is in target promoter region using GeneDatabase.

        Args:
            position_nt: List of nucleotide positions

        Returns:
            True if variant is within the target promoter region
        """
        promoter_region = GeneDatabase.get_promoter_region(self.gene_id)
        return Helper.is_mutation_within_range(position_nt, promoter_region)

    def _is_loss_of_function(self) -> bool:
        """Check if mutation represents loss-of-function per rule 2.2.1.1.

        Loss-of-function mutations contain: del, ins, fs, delins, _ or ends with *
        """
        lof_indicators = ["del", "ins", "fs", "delins", "_"]
        mutation_type = [self.nucleotide_change, self.protein_change]

        # Check both nucleotide and protein changes for LOF indicators
        for mutation in mutation_type:
            if any(indicator in mutation for indicator in lof_indicators) or mutation.endswith("*"):
                return True
        return False

    def _is_in_orf(self) -> bool:
        """Check if mutation is in the ORF.

        Used for rule 2.2.1.1 (LOF expert rule).
        """
        position_nt = Helper.get_position(self.nucleotide_change)
        if position_nt == [None]:
            return False
        return any(p >= 1 for p in position_nt)

    def _is_upstream_30_bp(self) -> bool:
        """Check if mutation is in ORF or within first 30 nucleotides upstream of start codon.

        Used for rule 2.2.1.1 (LOF expert rule).
        """
        position_nt = Helper.get_position(self.nucleotide_change)
        return any(int(position) > -30 for position in position_nt)

    def _is_deletion_in_orf(self) -> bool:
        """Check if a variant represents a deletion.

        Returns:
            True if the nucleotide_change contains 'del'
        """
        return "del" in self.nucleotide_change and self._is_in_orf()