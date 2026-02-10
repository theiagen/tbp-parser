from pydantic import BaseModel, computed_field

class GeneCoverage(BaseModel):
    """
      A class to represent coverage data for a single gene (BedRecord).
    """
    locus_tag: str
    gene_name: str
    coords: tuple[int, int]
    breadth_of_coverage: float
    average_depth: float

    def __str__(self):
        return f"GeneCoverage([{self.gene_name}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def __repr__(self):
        return f"GeneCoverage([{self.gene_name}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

class LocusCoverage(BaseModel):
    """
    A class to represent coverage data for a single locus tag, aggregating multiple regions if necessary.
    """
    locus_tag: str
    gene_names: list[str]
    coords: list[tuple[int, int]] # can be a list of coordinates if aggregating multiple regions
    breadth_of_coverage: float
    average_depth: float

    def __str__(self):
        return f"LocusCoverage([{self.gene_names}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def __repr__(self):
        return f"LocusCoverage([{self.gene_names}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def contains_position(self, position: int) -> bool:
        """Check if a specific position exists within locus coverage.

        Args:
            position (int): The position to check.
        Returns:
            bool: True if the position is within any of the coordinate ranges, False otherwise.
        """
        for start, end in self.coords:
            if start <= position <= end:
                return True
        return False

    def has_breadth_below(self, threshold: float) -> bool:
        """Check if the locus coverage has low breadth of coverage.

        Args:
            threshold (float): The breadth of coverage threshold.
        Returns:
            bool: True if the breadth of coverage is below the threshold, False otherwise.
        """
        return self.breadth_of_coverage < threshold * 100