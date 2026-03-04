from typing import Optional
from pydantic import BaseModel, Field, field_validator
from abc import ABC
from Variant import Variant

class BaseCoverage(BaseModel, ABC):
    """
    A base class to represent coverage data for a genomic region.
    This class can be extended to represent coverage for specific types of regions (e.g., locus/target/ERR).
    """
    coords: list[tuple[int, int]]
    breadth_of_coverage: float
    average_depth: float
    valid_deletions: list[Variant] = Field(default_factory=list) # list of Variant objects with valid deletions that fall within the coverage region

    def _coverage_str(self) -> str:
        return f"[BC:{(self.breadth_of_coverage * 100):.3f}%, AD:{self.average_depth:.3f}]"

    def __str__(self):
        return f"{self.__class__.__name__}({self._coverage_str()})"

    def __repr__(self):
        return self.__str__()

    @field_validator('breadth_of_coverage')
    @classmethod
    def check_breadth(cls, v):
        if v > 1.0:
            raise ValueError(f"breadth_of_coverage cannot exceed 1.0, got {v}")
        return v

    def has_breadth_below(self, threshold: float) -> bool:
        return self.breadth_of_coverage < threshold

    def contains_valid_deletion(self, variant: Variant) -> bool:
        """Check if a given Variant with a deletion is in the list of valid deletions for this coverage region."""
        return any(v == variant for v in self.valid_deletions)

    def contains_position(self, position: int) -> bool:
        """Check if a specific position exists within this coverage region."""
        return any(start <= position <= end for start, end in self.coords)

    def overlaps_range(self, start: int, end: int) -> bool:
        """Check if a genomic range [start, end] overlaps with this coverage region."""
        return any(cs <= end and start <= ce for cs, ce in self.coords)

class ERRCoverage(BaseCoverage):
    """
    A class to represent coverage data for the essential reportable range (ERR) of a single target/locus.
    """
    coords: list[tuple[int, int]]
    model_config = {"extra": "ignore"}

class TargetCoverage(BaseCoverage):
    """
      A class to represent coverage data for a single target genomic region (BedRecord).
      Based on unique gene_name's in BedRecord. Multiple TargetCoverage objects may exist for a single locus_tag/gene
      if there are multiple target regions (BedRecords) within that locus_tag/gene.
    """
    locus_tag: str
    gene_name: str
    coords: list[tuple[int, int]]
    err_coverage: Optional[ERRCoverage] = None

    def model_post_init(self, __context) -> None:
        if self.err_coverage is not None:
            if not all(self.contains_position(s) and self.contains_position(e) for s, e in self.err_coverage.coords):
                raise ValueError(f"ERR coords {self.err_coverage.coords} fall outside target coords {self.coords}")

class LocusCoverage(BaseCoverage):
    """
    A class to represent coverage data for a single locus tag, aggregating multiple target regions if necessary.
    """
    locus_tag: str
    gene_names: list[str]
    coords: list[tuple[int, int]] # can be a list of coordinates if aggregating multiple target regions
    err_coverage: Optional[ERRCoverage] = None

    # Post-init validation to ensure that if ERR coverage is provided, its coordinates fall within the locus coverage coordinates
    def model_post_init(self, __context) -> None:
        if self.err_coverage is not None:
            if not all(self.contains_position(s) and self.contains_position(e) for s, e in self.err_coverage.coords):
                raise ValueError(f"ERR coords {self.err_coverage.coords} fall outside locus coords {self.coords}")
