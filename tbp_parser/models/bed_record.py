from typing import Optional, Dict, List

class BedRecord:
    """A class representing a record or entry from a BED file."""

    def __init__(self, chrom: str, start: int, end: int, name: str, locus_tag: str, drug: list[str]) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.locus_tag = locus_tag
        self.name = name
        self.drug = drug
        self.reads_by_position: Dict[int, List[str]] = {} # to be populated in Coverage class (0-based)
        self.breadth_of_coverage: Optional[float] = None
        self.average_depth: Optional[float] = None

    def __str__(self):
        return f"BedRecord([{self.name}][{self.locus_tag}]{self.coords})"

    def __len__(self) -> int:
        return self.end - self.start + 1  # assuming 1-based indexing

    @property
    def length(self) -> int:
        return self.end - self.start + 1  # assuming 1-based indexing

    @property
    def coords(self) -> list[int]:
        return [self.start, self.end]

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
            name=cols[4],
            drug=cols[5].strip().split(',')
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

    def overlapping_coords(self, other: 'BedRecord') -> list[int]:
        """Get the overlapping coordinates between this BedRecord and another BedRecord.

        Args:
            other (BedRecord): Another BedRecord to get overlapping coordinates with.
        Returns:
            list[int]: A list containing the start and end of the overlapping region, or an empty list if there is no overlap.
        """
        if not self.overlaps_with(other):
            return []
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        return [overlap_start, overlap_end]

    def get_non_overlapping_reads(self, overlapping_bed_records: list['BedRecord']) -> set[str]:
        """Get a list of unique reads that cover the specified range within this BedRecord.
        Needed for when a bed_record overlaps with one or more other bed_records to avoid double-counting reads.

        Args:
            overlapping_bed_records (list['BedRecord']): A list of BedRecord instances that overlap with this BedRecord.
        Returns:
            set[str]: A set of unique read names covering the specified range.
        """
        # collect all positions from this BedRecord
        all_positions = set(range(self.start, self.end + 1))

        # remove positions that overlap with other bed_records
        for other in overlapping_bed_records:
            overlap_coords = self.overlapping_coords(other)
            if overlap_coords:
                overlap_start, overlap_end = overlap_coords
                # remove overlapping positions from the set
                all_positions -= set(range(overlap_start, overlap_end + 1))

        # get reads from non-overlapping positions only in this BedRecord
        unique_reads = set()
        for pos in all_positions:
            if pos in self.reads_by_position:
                unique_reads.update(self.reads_by_position[pos])
        return unique_reads