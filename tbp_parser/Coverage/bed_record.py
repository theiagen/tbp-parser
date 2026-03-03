from typing import Any, Dict, List
from pydantic import BaseModel, Field
from Utilities import Configuration

class BedRecord(BaseModel):
    """A class representing a record or entry from a BED file."""
    chrom: str
    start: int
    end: int
    locus_tag: str
    gene_name: str

    # Derived fields (computed during init, excluded from serialization)
    length: int = Field(default=0, exclude=True)
    coords: tuple[int, int] = Field(default=(0, 0), exclude=True)

    # To be populated in Coverage class after parsing the BAM file, excluded from serialization
    reads_by_position: Dict[int, List[str]] = Field(default_factory=dict, exclude=True) # (1-based)

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        # Calculate length and coords
        self.length = self.end - self.start + 1  # assuming 1-based indexing
        self.coords = (self.start, self.end)

        # Normalize gene_name
        Configuration.get_instance().normalize_field_values(self)

    def __str__(self):
        return f"BedRecord([{self.gene_name}][{self.locus_tag}]{self.coords})"

    def __repr__(self):
        return f"BedRecord([{self.gene_name}][{self.locus_tag}]{self.coords})"

    def __eq__(self, other):
        """Define equality based on some attributes"""
        if not isinstance(other, BedRecord):
            return False
        return (
            self.chrom == other.chrom and
            self.start == other.start and
            self.end == other.end and
            self.locus_tag == other.locus_tag and
            self.gene_name == other.gene_name
        )

    def __hash__(self):
        """Make the object hashable using some attributes"""
        return hash((
            self.chrom,
            self.start,
            self.end,
            self.locus_tag,
            self.gene_name
        ))

    @classmethod
    def from_bed_line(cls, bed_line: str) -> 'BedRecord':
        """Create a BedRecord instance from a tab separated line in a BED file.

        Args:
            bed_line (str): A line from a BED file.
        Returns:
            BedRecord: An instance of BedRecord.
        """
        cols = bed_line.strip().split('\t')
        return cls(
            chrom=cols[0],
            start=int(cols[1]),
            end=int(cols[2]),
            locus_tag=cols[3],
            gene_name=cols[4],
        )

    def overlaps_with(self, other: 'BedRecord') -> bool:
        """Check if this BedRecord overlaps with another BedRecord.

        Args:
            other (BedRecord): Another BedRecord to check overlap with.
        Returns:
            bool: True if there is an overlap, False otherwise.
        """
        overlap = (min(self.end, other.end) - max(self.start, other.start)) >= 0
        return overlap

    def overlapping_coords(self, other: 'BedRecord') -> tuple[int, int]:
        """Get the overlapping coordinates between this BedRecord and another BedRecord.

        Args:
            other (BedRecord): Another BedRecord to get overlapping coordinates with.
        Returns:
            tuple[int, int]: A tuple containing the start and end of the overlapping region.
        """
        # exception here because this should only be called if an overlap exists
        # can change this behavior if needed later
        if not self.overlaps_with(other):
            raise Exception("No overlap exists between BedRecords")
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        return (overlap_start, overlap_end)

    def get_non_overlapping_coords(self, others: list['BedRecord']) -> list[tuple[int, int]]:
        """Get the non-overlapping coordinates between this BedRecord and another BedRecord.
        Args:
            others (list['BedRecord']): A list of BedRecord instances that overlap with this BedRecord.
        Returns:
            list[tuple[int, int]]: A list of tuples containing the start and end of non-overlapping regions.
                                   Empty list when self is completely contained within other.
                                   1 tuple when they partially overlap on one side.
                                   2 tuples when other is completely contained within self.
        """
        non_overlapping_positions = self._get_non_overlapping_positions(others)
        if not non_overlapping_positions:
            return []

        # convert set of ordered positions back to a list of tuples
        coords = []
        sorted_positions = sorted(non_overlapping_positions)
        start = sorted_positions[0]
        prev = start

        # if a gap is detected, that indicates the end of a non-overlapping region
        for pos in sorted_positions[1:]:
            if pos != prev + 1:
                coords.append((start, prev))
                start = pos
            prev = pos

        # add the last region
        coords.append((start, prev))

        # sanity check - should be at most 2 non-overlapping regions
        if len(coords) > 2:
            raise ValueError(f"Expected at most 2 non-overlapping regions, got {len(coords)}: {coords}")
        return coords

    def get_non_overlapping_reads(self, others: list['BedRecord']) -> set[str]:
        """Get a list of unique reads that cover the specified range within this BedRecord.
        Needed for when a bed_record overlaps with one or more other bed_records to avoid double-counting reads.

        Args:
            others (list['BedRecord']): A list of BedRecord instances that overlap with this BedRecord.
        Returns:
            set[str]: A set of unique read names covering the specified range.
        """
        non_overlapping_positions = self._get_non_overlapping_positions(others)

        # get reads from non-overlapping positions only in this BedRecord
        unique_reads = set()
        for pos in non_overlapping_positions:
            if pos in self.reads_by_position:
                unique_reads.update(self.reads_by_position[pos])
        return unique_reads

    def _get_non_overlapping_positions(self, others: list['BedRecord']) -> set[int]:
        """Get the non-overlapping positions between this BedRecord and another BedRecord.

        Args:
            others (list['BedRecord']): A list of BedRecord instances that overlap with this BedRecord.
        Returns:
            set[int]: A set of non-overlapping positions.
        """
        # collect all positions from this BedRecord
        non_overlapping_positions = set(range(self.start, self.end + 1))

        # remove positions that overlap with `other` bed_records
        for other in others:
            if self.overlaps_with(other):
                overlap_start, overlap_end = self.overlapping_coords(other)
                # remove overlapping positions from the set
                non_overlapping_positions -= set(range(overlap_start, overlap_end + 1))

        return non_overlapping_positions