import logging
import pysam
from collections import defaultdict
from itertools import combinations
from typing import Dict, Optional, List
from collections import defaultdict

from utils.config import Configuration
from models.bed_record import BedRecord
from models.coverage_data import GeneCoverage, LocusCoverage

logger = logging.getLogger(__name__)

class CoverageCalculator:
    """A class used to generage coverage statistics for BED records."""

    def __init__(self, config: Configuration) -> None:
        self.config = config

    @staticmethod
    def _fetch_reads_by_position(
        bed_record: BedRecord,
        input_bam: pysam.AlignmentFile,
        whitelisted_reads: Optional[set[str]] = None
    ) -> Dict[int, List[str]]:
        """
        Calculate coverage for a single BedRecord.

        Args:
            bed_record: The BedRecord to calculate coverage for
            whitelisted_reads: Optional set of read names to filter by (for overlap handling)
        """
        reads_by_position = {}
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
                #     logger.debug(f"{bed_record} | Position: {rec.pos + 1} | Total Reads: {len(rec.get_query_names())} | Non-overlapping Reads: {len(read_list)} | Reads Removed: {reads_removed}")

            reads_by_position[rec.pos + 1] = read_list # convert to 1-based indexing
        return reads_by_position


    def populate_reads_by_position(self, bed_records: List[BedRecord]) -> List[BedRecord]:
        """
        Populates the reads_by_position attribute for all BedRecords using the config BAM file.
        Args:
            bed_records (List[BedRecord]): The list of BedRecords to populate reads_by_position
        """
        logger.debug(f"Populating reads_by_position for {len(bed_records)} BedRecords")

        with pysam.AlignmentFile(self.config.input_bam, "rb") as input_bam:
            for bed_record in bed_records:
                bed_record.reads_by_position = self._fetch_reads_by_position(bed_record, input_bam)
                logger.debug(f"Populated reads_by_position for {bed_record} with {len(bed_record.reads_by_position)} positions")
        return bed_records


    def _calculate_breadth_of_coverage(self, reads_by_position: dict[int, list[str]], length: int) -> float:
        # Calculates breadth of coverage for a given BedRecord.
        positions_above_min_depth = sum(1 for reads in reads_by_position.values() if len(reads) >= self.config.MIN_DEPTH)
        return (positions_above_min_depth / length) * 100


    def _calculate_average_depth(self, reads_by_position: dict[int, list[str]], length: int) -> float:
        # Calculates average depth of coverage for a given BedRecord.
        depth_sum = sum(len(reads) for reads in reads_by_position.values())
        return depth_sum / length


    def resolve_overlapping_regions(self, bed_records: List[BedRecord]) -> List[BedRecord]:
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
        logger.info(f"Found {len(overlap_map)} overlapping records, reassigning BedRecord.reads_by_position with non-overlapping regions only")

        with pysam.AlignmentFile(self.config.input_bam, "rb") as input_bam:
            for bed_record in overlap_map.keys():
                overlapping_bed_records = overlap_map[bed_record]
                logger.debug(f"Resolving overlapping region for {bed_record} which overlaps with {len(overlapping_bed_records)} other BedRecord(s)")

                # get names of all reads in non-overlapping regions for this bed_record
                whitelisted_reads = bed_record.get_non_overlapping_reads(overlapping_bed_records)

                # just for logging/debugging
                non_overlapping_coords = bed_record.get_non_overlapping_coords(overlapping_bed_records)
                logger.debug(f"{bed_record} | Non-overlapping coords: {non_overlapping_coords}")

                # reassigning reads_by_position for this bedrecord only considering the whitelisted reads
                bed_record.reads_by_position = self._fetch_reads_by_position(bed_record, input_bam, whitelisted_reads)
        return bed_records


    def generate_gene_coverage_list(self, bed_records: List[BedRecord]) -> List[GeneCoverage]:
        """
        Calculates gene-level coverage statistics for the provided BedRecords.
        One GeneCoverage object is created per BedRecord.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to calculate gene coverage for.
        Returns:
            List[GeneCoverage]: A list of GeneCoverage objects representing coverage statistics for each gene.
        """
        gene_coverage_list = []
        for bed_record in bed_records:
            total_length = bed_record.length
            breadth_of_coverage = self._calculate_breadth_of_coverage(bed_record.reads_by_position, total_length)
            average_depth = self._calculate_average_depth(bed_record.reads_by_position, total_length)

            gene_coverage = GeneCoverage(
                locus_tag=bed_record.locus_tag,
                gene_name=bed_record.gene_name,
                coords=bed_record.coords,
                breadth_of_coverage=breadth_of_coverage,
                average_depth=average_depth,
            )
            gene_coverage_list.append(gene_coverage)
        return gene_coverage_list


    def generate_locus_coverage_list(self, bed_records: List[BedRecord]) -> List[LocusCoverage]:
        """
        Calculates locus-level coverage statistics for the provided BedRecords.
        Aggregates multiple BedRecords with the same locus tag into a single LocusCoverage object.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to calculate locus coverage for.
        Returns:
            List[LocusCoverage]: A list of LocusCoverage objects representing coverage statistics for each locus tag.
        """

        locus_groups = defaultdict(list)
        for bed_record in bed_records:
            locus_groups[bed_record.locus_tag].append(bed_record)

        # at this point, overlaps should be resolved, so we can just aggregate reads_by_position
        locus_coverage_list = []
        for locus_tag, bed_records in locus_groups.items():
            all_reads_by_position = defaultdict(list)
            for bed_record in bed_records:
                for pos, reads in bed_record.reads_by_position.items():
                    all_reads_by_position[pos].extend(reads)

            total_length = len(all_reads_by_position)
            breadth_of_coverage = self._calculate_breadth_of_coverage(all_reads_by_position, total_length)
            average_depth = self._calculate_average_depth(all_reads_by_position, total_length)

            locus_coverage = LocusCoverage(
                locus_tag=locus_tag,
                gene_names=[r.gene_name for r in bed_records],
                coords=[r.coords for r in bed_records],
                breadth_of_coverage=breadth_of_coverage,
                average_depth=average_depth,
            )
            locus_coverage_list.append(locus_coverage)
        return locus_coverage_list