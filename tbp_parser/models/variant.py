import copy
import globals as globals_
from Row import Row
from typing import Optional
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
        pos: int,
        depth: int,
        freq: float,
        gene_id: str,
        gene_name: str,
        type: str,
        nucleotide_change: str,
        confidence: str,
        drug: str,
        source: str,
        comment: str,
        protein_change: str = "",
    ) -> None:
        """Initializes the Variant class

        Args:
            variant_data (dict): a dictionary of attributes about the variant.
        """

        # initialize the variant attributes from variant_data dictionary. purposely error if key is missing
        self.pos: int = pos
        self.depth: int = depth
        self.freq: float = freq
        self.gene_id: str = gene_id
        self.gene_name: str = gene_name
        self.type: str = type
        self.nucleotide_change: str = nucleotide_change
        # annotation-related attributes
        self.confidence = confidence
        self.drug = drug
        self.source = source
        self.comment = comment
        # private attribute to access via property
        self._protein_change: str = protein_change

        # QC and interpretation-related attributes to be populated later
        self.rationale: Optional[str] = None
        self.looker_interpretation: Optional[str] = None
        self.mdl_interpretation: Optional[str] = None
        self.warning: set[str] = set()

    @property
    def protein_change(self) -> str:
        """Helper method to get the `change` attribute based on variant type"""
        if self.type == "synonymous_variant" or not self._protein_change:
            return self.nucleotide_change
        return self._protein_change

    @property
    def read_support(self) -> float:
        return self.freq * self.depth


    @classmethod
    def from_thin_air(cls, gene_id: str, gene_name: str, drug: str) -> 'Variant':
        """Creates a Variant object with minimal information for unreported/WT variants.

        Args:
            gene_id (str): The gene ID of the variant.
            gene_name (str): The gene name of the variant.
            drug (str): The drug associated with the variant.
        Returns:
            Variant: A Variant object with minimal information.
        """
        return cls(
            pos=-1,
            depth=-1,
            freq=-1.0,
            gene_id=gene_id,
            gene_name=gene_name,
            type="NA",
            nucleotide_change="NA",
            protein_change="NA",
            confidence="NA",
            drug=drug,
            source="",
            comment="",
        )

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
          "Assoc w R - Interim": 4,
          "Uncertain significance": 3,
          "Not assoc w R - Interim": 2,
          "Not assoc w R": 1,
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

