from typing import Optional
from utils import GeneDatabase, Helper

import logging

logger = logging.getLogger(__name__)

class Variant:
    """This class represents a single expanded variant reported by TBProfiler.

    This class assumes that the variant has already been expanded such that
    each Variant instance represents a single gene-drug association. These
    associations can be derived from the annotation and consequence field(s)
    in the tbp-profiler output JSON.

    """
    def __init__(
        self,
        sample_id: str,
        pos: int,
        depth: int,
        freq: float,
        gene_id: str,
        gene_name: str,
        type: str,
        nucleotide_change: str,
        protein_change: str,
        confidence: str,
        drug: str,
        source: str,
        comment: str,
    ) -> None:
        """Initializes the Variant class

        Args:
            variant_data (dict): a dictionary of attributes about the variant.
        """

        # initialize the variant attributes from input JSON dictionary. purposely error if key is missing
        self.sample_id: str = sample_id
        self.pos: int = pos
        self.depth: int = depth
        self.freq: float = freq
        self.read_support: float = self.freq * self.depth
        self.gene_id: str = gene_id
        self.gene_name: str = gene_name
        self.type: str = type
        self.nucleotide_change: str = nucleotide_change
        self.protein_change: str = protein_change

        # annotation-related attributes
        self.confidence: str = confidence
        self.drug: str = drug
        self.source: str = source
        self.comment: str = comment

        # VariantInterpreter will assign these attributes later
        self.rationale: Optional[str] = None
        self.looker_interpretation: Optional[str] = None
        self.mdl_interpretation: Optional[str] = None

        # VariantQC will assign these attributes later
        self.fails_qc: Optional[bool] = False
        self.warning: set[str] = set()

        # derive these attributes from various functions
        self.gene_tier: str = GeneDatabase.get_tier(self.gene_name) if self.gene_name in GeneDatabase._GENE_NAME_TO_LOCUS else "NA"
        self.absolute_start: int = Helper.get_mutation_genomic_positions(self.pos, self.nucleotide_change)[0]
        self.absolute_end: int = Helper.get_mutation_genomic_positions(self.pos, self.nucleotide_change)[1]

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

    def is_better_annotation_than(self, other: 'Variant') -> bool:
        annotation_rank = {
          "Assoc w R": 5,
          "Assoc w R - interim": 4,
          "Assoc w R - Interim": 4,
          "Uncertain significance": 3,
          "Not assoc w R - Interim": 2, # should these be flipped?
          "Not assoc w R": 1, # should these be flipped?
          "No WHO annotation": 0,
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