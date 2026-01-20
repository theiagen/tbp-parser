import logging
import pysam
from collections import defaultdict
from itertools import combinations
from typing import Optional, List

from utils.config import Configuration
from models.bed_record import BedRecord

logger = logging.getLogger(__name__)

class CoverageCalculator:
    """A class used to generage coverage statistics for BED records."""

    def __init__(self, config: Configuration) -> None:
        self.config = config


    def generate_cov_attr(self, bed_records: List[BedRecord]) -> List[BedRecord]:
        """
        Populates the coverage attributes for all BedRecords using the config BAM file.

        Args:
            bed_records (List[BedRecord]): The list of BedRecords to populate coverage attributes for.
        Returns:
            List[BedRecord]: The list of BedRecords with populated coverage attributes.
        """
        logger.debug(f"Populating coverage attributes for {len(bed_records)} BedRecords")

        with pysam.AlignmentFile(self.config.input_bam, "rb") as input_bam:
            for bed_record in bed_records:
                self._calculate_coverage(bed_record, input_bam)
                logger.debug(f"{bed_record}: Breadth of Coverage: {bed_record.breadth_of_coverage}, Average Depth: {bed_record.average_depth}")
        return bed_records


    def _calculate_coverage(
        self,
        bed_record: BedRecord,
        input_bam: pysam.AlignmentFile,
        whitelisted_reads: Optional[set[str]] = None) -> None:
        """
        Calculate coverage for a single BedRecord.

        Args:
            bed_record: The BedRecord to calculate coverage for
            whitelisted_reads: Optional set of read names to filter by (for overlap handling)
        """
        bam_chrom = input_bam.get_reference_name(0) # just in case theres a chance this chromosome name is different from the bed file
        for rec in input_bam.pileup(
            contig=bam_chrom,
            start=bed_record.start - 1, # start uses 0-based indexing in pysam
            end=bed_record.end, # end is non-inclusive in pysam
            stepper="nofilter", # include all reads mimics samtools depth as closely as possible (not exact)
            truncate=True,
        ):
            read_list = rec.get_query_names()

            if whitelisted_reads:
                read_list = [read for read in read_list if read in whitelisted_reads]
                # reads_removed = len(rec.get_query_names()) - len(read_list)
                # if reads_removed > 0:
                #     self.logger.debug(f"{bed_record}Position: {rec.pos + 1} | Total Reads: {len(rec.get_query_names())} | Non-overlapping Reads: {len(read_list)} | Reads Removed: {reads_removed}")

            bed_record.reads_by_position[rec.pos + 1] = read_list # convert to 1-based indexing

        # Calculates coverage statistics for a given BedRecord.
        total_positions = bed_record.length
        positions_above_min_depth = sum(1 for reads in bed_record.reads_by_position.values() if len(reads) >= self.config.MIN_DEPTH)
        depth_sum = sum(len(reads) for reads in bed_record.reads_by_position.values())

        bed_record.breadth_of_coverage = (positions_above_min_depth / total_positions) * 100
        bed_record.average_depth = depth_sum / total_positions


    def update_overlapping_bed_records(self, bed_records: List[BedRecord]) -> List[BedRecord]:
        """
        Update coverage for all overlapping regions in the provided list of BedRecords.
        For records with overlaps, recalculates coverage using only reads
        in non-overlapping regions to avoid double-counting.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to check for overlaps and
                update coverage statistics.
        Returns:
            List[BedRecord]: The list of BedRecords with updated coverage statistics.
        """

        logger.info("Checking for overlapping regions in the BED file")

        overlap_map = defaultdict(list)
        # build overlap map: record -> list of records it overlaps with
        for bed_record1, bed_record2 in combinations(bed_records, 2):
            if bed_record1.overlaps_with(bed_record2):
                overlap_map[bed_record1].append(bed_record2)
                overlap_map[bed_record2].append(bed_record1)

        if not overlap_map:
            logger.info("No overlapping regions found in the BED file")
            return bed_records

        # found overlaps - need to recalculate
        logger.info(f"Found {len(overlap_map)} overlapping records, recalculating coverage")

        with pysam.AlignmentFile(self.config.input_bam, "rb") as input_bam:
            for bed_record in overlap_map.keys():
                overlapping_bed_records = overlap_map[bed_record]
                logger.debug(f"Found overlapping region for {bed_record} which overlaps with {len(overlapping_bed_records)} other BedRecord(s)")

                # get names of all reads in non-overlapping regions for this bed_record
                whitelisted_reads = bed_record.get_non_overlapping_reads(overlapping_bed_records)

                # store old values just for logging
                old_breadth = bed_record.breadth_of_coverage
                old_depth = bed_record.average_depth

                # recalculate coverage attributes of this bedrecord only considering the whitelisted reads
                self._calculate_coverage(bed_record, input_bam, whitelisted_reads)

                logger.debug(f"BEFORE: {bed_record}: Breadth of Coverage: {old_breadth}, Average Depth: {old_depth}")
                logger.debug(f"AFTER: {bed_record}: Breadth of Coverage: {bed_record.breadth_of_coverage}, Average Depth: {bed_record.average_depth}")
        return bed_records


    def get_breadth_of_coverage_map(self, bed_records: List[BedRecord]) -> dict[str, float]:
        """
        Generates a coverage map for the provided BedRecords.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to generate the coverage map for.
        Returns:
            dict[str, float]: A dictionary mapping BedRecord names to their breadth of coverage.
        """
        coverage_map = {}
        for bed_record in bed_records:
            coverage_map[bed_record.name] = bed_record.breadth_of_coverage
        return coverage_map


    def get_depth_of_coverage_map(self, bed_records: List[BedRecord]) -> dict[str, float]:
        """
        Generates a depth of coverage map for the provided BedRecords.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to generate the depth of coverage map for.
        Returns:
            dict[str, float]: A dictionary mapping BedRecord names to their average depth of coverage.
        """
        coverage_map = {}
        for bed_record in bed_records:
            coverage_map[bed_record.name] = bed_record.average_depth
        return coverage_map


    def get_low_coverage_list(self, bed_records: List[BedRecord]) -> List[str]:
        """
        Identifies BedRecords with breadth of coverage below the minimum threshold.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to check for low coverage.
        Returns:
            List[str]: A list of BedRecord names with breadth of coverage below the minimum threshold.
        """
        low_coverage_records = []
        for bed_record in bed_records:
            if bed_record.breadth_of_coverage < (self.config.MIN_PERCENT_COVERAGE * 100):
                low_coverage_records.append(bed_record.name)
        return low_coverage_records