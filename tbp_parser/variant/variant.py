from pydantic import BaseModel, Field
from typing import Optional, Any
from utils import GeneDatabase, Helper
import logging

logger = logging.getLogger(__name__)

class Variant(BaseModel):
    """This class represents a single genetic variant that can be found in the Laboratorian report.

    Core attributes include fields expanded from the TBProfiler output JSON.

    This class assumes that the variant has already been expanded such that
    each Variant instance represents a single gene-drug association. These
    associations can be derived from the annotation and consequence field(s)
    in the tbp-profiler output JSON.

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
    fails_qc: Optional[bool] = Field(default=None, exclude=True)
    warning: set[str] = Field(default_factory=set, exclude=True)

    model_config = {"extra": "ignore"}

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        # Calculate read_support if not provided
        if self.read_support is None:
            self.read_support = self.freq * self.depth

        # Derive computed attributes
        self.gene_tier = GeneDatabase.get_tier(self.gene_name)
        self.absolute_start, self.absolute_end = Helper.get_mutation_genomic_positions(self.pos, self.nucleotide_change)

        # Normalize field values based on predefined rules
        Helper.normalize_field_values(self)

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
          "Not assoc w R - Interim": 2, # should these be flipped?
          "Not assoc w R": 1, # should these be flipped?
          "Not found in WHO catalogue": 0, # this might be redundant with "No WHO annotation"
          "No WHO annotation": -1, # given to synthetic variants
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