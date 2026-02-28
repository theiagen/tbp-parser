from typing import Optional
from pydantic import BaseModel, Field, field_validator
from abc import ABC, abstractmethod
from variant import Variant

class BaseCoverage(BaseModel, ABC):
    """
    A base class to represent coverage data for a genomic region.
    This class can be extended to represent coverage for specific types of regions (e.g., locus/target/ERR).
    """
    breadth_of_coverage: float
    average_depth: float
    valid_deletions: list[Variant] = Field(default_factory=list) # list of Variant objects with valid deletions that fall within the coverage region

    def _coverage_str(self) -> str:
        return f"[BC:{(self.breadth_of_coverage * 100):.3f}%, AD:{self.average_depth:.3f}]"

    def __str__(self):
        return f"{self.__class__.__name__}({self._coverage_str()})"

    def __repr__(self):
        return self.__str__()

    @abstractmethod
    def contains_position(self, position: int) -> bool: ...

    @field_validator('breadth_of_coverage')
    @classmethod
    def check_breadth(cls, v):
        if v > 1.0:
            raise ValueError(f"breadth_of_coverage cannot exceed 1.0, got {v}")
        return v

    def has_breadth_below(self, threshold: float) -> bool:
        return self.breadth_of_coverage < threshold

class ERRCoverage(BaseCoverage):
    """
    A class to represent coverage data for the essential reportable range (ERR) of a single target/locus.
    """
    coords: tuple[int, int] | list[tuple[int, int]]
    model_config = {"extra": "ignore"}

    def contains_position(self, position: int) -> bool:
        """Check if a specific position exists within the ERR coverage.
        Args:
            position (int): The position to check.
          Returns:
              bool: True if the position is within any of the coordinate ranges, False otherwise.
        """
        normalized = self.coords if isinstance(self.coords, list) else [self.coords]
        return any(start <= position <= end for start, end in normalized)

class TargetCoverage(BaseCoverage):
    """
      A class to represent coverage data for a single target genomic region (BedRecord).
      Based on unique gene_name's in BedRecord. Multiple TargetCoverage objects may exist for a single locus_tag/gene
      if there are multiple target regions (BedRecords) within that locus_tag/gene.
    """
    locus_tag: str
    gene_name: str
    coords: tuple[int, int]
    err_coverage: Optional[ERRCoverage] = None

    def contains_position(self, position: int) -> bool:
        """Check if a specific position exists within target coverage.

        Args:
            position (int): The position to check.
        Returns:
            bool: True if the position is within any of the coordinate ranges, False otherwise.
        """
        return self.coords[0] <= position <= self.coords[1]

class LocusCoverage(BaseCoverage):
    """
    A class to represent coverage data for a single locus tag, aggregating multiple target regions if necessary.
    """
    locus_tag: str
    gene_names: list[str]
    coords: list[tuple[int, int]] # can be a list of coordinates if aggregating multiple target regions
    err_coverage: Optional[ERRCoverage] = None

    def contains_position(self, position: int) -> bool:
        """Check if a specific position exists within locus coverage.

        Args:
            position (int): The position to check.
        Returns:
            bool: True if the position is within any of the coordinate ranges, False otherwise.
        """
        return any(start <= position <= end for start, end in self.coords)