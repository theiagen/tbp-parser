
class GeneCoverage:
    """
      A class to represent coverage data for a single gene (BedRecord).
    """
    def __init__(
        self,
        locus_tag: str,
        gene_name: str,
        coords: tuple[int, int],
        breadth_of_coverage: float,
        average_depth: float,
    ) -> None:

        self.locus_tag = locus_tag
        self.gene_name = gene_name
        self.coords = coords
        self.breadth_of_coverage = breadth_of_coverage
        self.average_depth = average_depth

    def __str__(self):
        return f"GeneCoverage([{self.gene_name}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def __repr__(self):
        return f"GeneCoverage([{self.gene_name}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def __len__(self) -> int:
        return self.length

    @property
    def length(self) -> int:
        return self.coords[1] - self.coords[0] # already 0-based from BedRecord


class LocusCoverage:
    """
    A class to represent coverage data for a single locus tag, aggregating multiple regions if necessary.
    """
    def __init__(
        self,
        locus_tag: str,
        gene_names: list[str],
        coords: list[tuple[int, int]], # can be a list of coordinates if aggregating multiple regions
        breadth_of_coverage: float,
        average_depth: float,
    ) -> None:

        self.locus_tag = locus_tag
        self.gene_names = gene_names
        self.coords = coords
        self.breadth_of_coverage = breadth_of_coverage
        self.average_depth = average_depth

    def __str__(self):
        return f"LocusCoverage([{self.gene_names}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def __repr__(self):
        return f"LocusCoverage([{self.gene_names}][{self.locus_tag}][{self.coords}] BC:{self.breadth_of_coverage} AD:{self.average_depth})"

    def __len__(self) -> int:
        return self.length

    @property
    def length(self) -> int:
        return sum(end - start for start, end in self.coords)

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