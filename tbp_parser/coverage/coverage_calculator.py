import logging
import pysam
from collections import defaultdict
from itertools import permutations
from typing import Dict, List, Tuple
from collections import defaultdict

from utils import Configuration
from coverage.bed_record import BedRecord
from coverage.coverage_data import TargetCoverage, LocusCoverage

logger = logging.getLogger(__name__)

class CoverageCalculator:
    """A class used to generage coverage statistics for BED records."""

    def __init__(self, config: Configuration) -> None:
        self.config = config


    def calculate(self, bed_records: List[BedRecord], err_records: List[BedRecord]) -> Tuple[Dict[str, LocusCoverage], Dict[str, TargetCoverage]]:
        """Main method to calculate coverage statistics for the provided BedRecords.

          Args:
              bed_records (List[BedRecord]): A list of BRR BedRecords to calculate coverage for.
              err_records (List[BedRecord]): A list of ERR BedRecords to calculate coverage for.
          Returns:
              Tuple[Dict[str, LocusCoverage], Dict[str, TargetCoverage]]: A tuple containing:
                  - Dictionary of LocusCoverage objects (one per locus tag)
                  - Dictionary of TargetCoverage objects (one per gene name)
        """

        # Fetch and populate all reads_by_position for all BedRecords first
        bed_records = self.populate_reads_by_position(bed_records)

        # Resolve overlaps if in tNGS mode to avoid double-counting reads across overlapping regions
        if self.config.TNGS:
            bed_records = self.resolve_overlapping_regions(bed_records)

        # Generate coverage maps
        locus_coverage_map, target_coverage_map = self.generate_coverage_maps(bed_records)

        # Call this function again with just err_records to populate reads, resolve overlaps, and generate coverage maps for ERR regions only
        # and set ERR coverage to each target and locus coverage record if applicable
        if err_records:
            logger.debug("Calculating ERR coverage and assigning to TargetCoverage and LocusCoverage records")
            err_locus_coverage_map, err_target_coverage_map = self.calculate(err_records, [])

            for gene_name, target_coverage in target_coverage_map.items():
                if gene_name in err_target_coverage_map:
                    setattr(target_coverage, "err_coverage", err_target_coverage_map[gene_name])

            for locus_tag, locus_coverage in locus_coverage_map.items():
                if locus_tag in err_locus_coverage_map:
                    setattr(locus_coverage, "err_coverage", err_locus_coverage_map[locus_tag])
        return locus_coverage_map, target_coverage_map


    def populate_reads_by_position(self, bed_records: List[BedRecord]) -> List[BedRecord]:
        """
        Populates the reads_by_position attribute for all BedRecords using the config BAM file.

        Notes and comparisons with samtools depth:
        See https://www.htslib.org/doc/samtools-depth.html#CAVEATS.
        See https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup.

        The following commands should produce identical read counts for each position:
        - samtools depth -a -J -r "Chromosome:100-200" input.bam
        - samtools mpileup -a --count-orphans --min-BQ 0 -r "Chromosome:100-200" input.bam
        - pysam AlignmentFile.pileup with stepper="nofilter", truncate=True, and the same region

        Args:
            bed_records: The list of BedRecords to calculate coverage for
        Returns:
            List[BedRecord]: The input list of BedRecords with the reads_by_position attribute populated for each record.
        """
        logger.debug(f"Populating reads_by_position for {len(bed_records)} BedRecords")

        with pysam.AlignmentFile(self.config.input_bam, "rb") as input_bam:
            for bed_record in bed_records:

                # Initialize with empty positions with no coverage. Pysam pileup only returns positions that have at least one read
                reads_by_position = {pos: [] for pos in range(bed_record.start, bed_record.end + 1)}
                bam_chrom = input_bam.get_reference_name(0) # just in case theres a chance this chromosome name is different from the bed file
                for rec in input_bam.pileup(
                    contig=bam_chrom,
                    start=bed_record.start - 1, # 0-based (inclusive)
                    end=bed_record.end, # 0-based (exclusive)
                    ignore_overlaps=False, # -s flag in samtools depth
                    ignore_orphans=False, # default behaviour in samtools depth
                    min_base_quality=0, # default in samtools depth
                    truncate=True,
                ):
                    read_list = rec.get_query_names()

                    # reads_removed = len(rec.get_query_names()) - len(read_list)
                    # if reads_removed > 0:
                    #     logger.debug(f"{bed_record} | Position: {rec.reference_pos + 1} | Total Reads: {len(rec.get_query_names())} | Non-overlapping Reads: {len(read_list)} | Reads Removed: {reads_removed}")

                    # convert to 1-based indexing.
                    # Ex) BedRecord (1-based closed interval) [1,10] -> pysam pileup [0,10)
                    reads_by_position[rec.reference_pos + 1] = read_list
                    setattr(bed_record, "reads_by_position", reads_by_position)
                logger.debug(f"Populated reads_by_position for {bed_record} across {len(bed_record.reads_by_position)} positions")
        return bed_records


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
        for bed_record1, bed_record2 in permutations(bed_records, 2):
            if (
                bed_record1.overlaps_with(bed_record2) and
                bed_record2 not in overlap_map[bed_record1]
            ):
                overlap_map[bed_record1].append(bed_record2)

        if not overlap_map:
            logger.info("No overlapping regions found in the BED file")
            return bed_records

        # found overlaps - need to recalculate
        logger.info(f"Found {len(overlap_map)} overlapping records, reassigning BedRecord.reads_by_position with non-overlapping regions only")

        # with pysam.AlignmentFile(self.config.input_bam, "rb") as input_bam:
        for bed_record in overlap_map.keys():
            overlapping_bed_records = overlap_map[bed_record]
            logger.debug(f"Resolving overlapping region for {bed_record} which overlaps with {len(overlapping_bed_records)} other BedRecord(s)")

            # get names of all reads in non-overlapping regions for this bed_record
            whitelisted_reads = bed_record.get_non_overlapping_reads(overlapping_bed_records)

            # just for logging/debugging
            non_overlapping_coords = bed_record.get_non_overlapping_coords(overlapping_bed_records)
            logger.debug(f"{bed_record} | Non-overlapping coords: {non_overlapping_coords}")

            # reassigning reads_by_position for this bedrecord only considering the whitelisted reads
            num_reads_kept = 0
            total_reads = 0
            for pos, read_list in bed_record.reads_by_position.items():
                total_reads += len(read_list)
                reads_kept = [read for read in read_list if read in whitelisted_reads]
                num_reads_kept += len(reads_kept)
                bed_record.reads_by_position[pos] = reads_kept
            logger.debug(f"Filtered and kept {num_reads_kept}/{total_reads} unique reads across all positions for {bed_record}")
        return bed_records


    def generate_coverage_maps(self, bed_records: List[BedRecord]) -> Tuple[Dict[str, LocusCoverage], Dict[str, TargetCoverage]]:
        """
        Calculates both target-level and locus-level coverage statistics for the provided BedRecords.
        For target-level coverage, one TargetCoverage object is created per BedRecord.
        For locus-level coverage (if applicable), multiple BedRecords with the same locus tag are aggregated into a single LocusCoverage object.

        Args:
            bed_records (List[BedRecord]): A list of BedRecords to calculate coverage for.
        Returns:
            Tuple[Dict[str, LocusCoverage], Dict[str, TargetCoverage]]: A tuple containing:
                - Dictionary of LocusCoverage objects (one per locus tag)
                - Dictionary of TargetCoverage objects (one per gene name)
        """
        target_coverage_map = {}
        locus_coverage_map = {}
        breadth_of_coverage = -1
        average_depth = -1

        # Group bed records by locus tag
        locus_groups = defaultdict(list)
        for bed_record in bed_records:
            locus_groups[bed_record.locus_tag].append(bed_record)

        # Process each locus group
        for locus_tag, bed_records in locus_groups.items():
            # Gene-level coverage (one per BedRecord)
            for bed_record in bed_records:
                logger.debug(f"Calculating breadth of coverage and average depth of coverage for gene: '{bed_record.gene_name}'")
                breadth_of_coverage = self._calculate_breadth_of_coverage(
                    bed_record.reads_by_position,
                )
                average_depth = self._calculate_average_depth(
                    bed_record.reads_by_position,
                )

                target_coverage_map[bed_record.gene_name] = TargetCoverage(
                    locus_tag=bed_record.locus_tag,
                    gene_name=bed_record.gene_name,
                    coords=bed_record.coords,
                    breadth_of_coverage=breadth_of_coverage,
                    average_depth=average_depth,
                )

            if len(bed_records) > 1:
                # Multiple BedRecords: aggregate reads_by_position
                logger.debug(f"Multiple BedRecords found for locus `{locus_tag}` ({[str(r.gene_name) for r in bed_records]}); Aggregating reads_by_position across {len(bed_records)} BedRecord(s)")

                all_reads_by_position = defaultdict(list)
                # no read is counted twice across multiple genes with the same locus tag since we resolved overlaps previously
                for bed_record in bed_records:
                    for pos, reads in bed_record.reads_by_position.items():
                        all_reads_by_position[pos].extend(reads)

                logger.debug(f"Calculating breadth of coverage and average depth of coverage for locus: `{locus_tag}` ({[str(r.gene_name) for r in bed_records]})")
                breadth_of_coverage = self._calculate_breadth_of_coverage(
                    all_reads_by_position,
                )
                average_depth = self._calculate_average_depth(
                    all_reads_by_position,
                )
            else:
                logger.debug(f"Single BedRecord found for locus `{locus_tag}` ({[str(r.gene_name) for r in bed_records]}), using previously calculated gene-level coverage statistics")

            # Locus-level coverage (one per locus tag, aggregating multiple BedRecords if necessary)
            # Otherwise, if only one BedRecord for this locus, the locus coverage will be the same as the gene coverage, just with a different object type
            locus_coverage = LocusCoverage(
                locus_tag=locus_tag,
                gene_names=[r.gene_name for r in bed_records],
                coords=[r.coords for r in bed_records],
                breadth_of_coverage=breadth_of_coverage,
                average_depth=average_depth,
            )
            locus_coverage_map[locus_tag] = locus_coverage
        return locus_coverage_map, target_coverage_map


    def _calculate_breadth_of_coverage(self, reads_by_position: dict[int, list[str]]) -> float:
        # Calculates breadth of coverage for a given BedRecord.
        positions_above_min_depth = sum(1 for reads in reads_by_position.values() if len(reads) >= self.config.MIN_DEPTH)
        try:
            breadth_of_coverage = positions_above_min_depth / len(reads_by_position)
            logger.debug(f"Total positions above MIN_DEPTH: {positions_above_min_depth}/{len(reads_by_position)} | Breadth of coverage: {(breadth_of_coverage * 100):.2f}%")
            return breadth_of_coverage
        except ZeroDivisionError:
            logger.debug("Length of BedRecord is 0, returning breadth of coverage as 0.0%")
            return 0.0


    def _calculate_average_depth(self, reads_by_position: dict[int, list[str]]) -> float:
        # Calculates average depth of coverage for a given BedRecord.
        depth_sum = sum(len(reads) for reads in reads_by_position.values())
        try:
            average_depth = depth_sum / len(reads_by_position)
            logger.debug(f"Total read depth across all positions: {depth_sum}/{len(reads_by_position)} | Average depth of coverage: {average_depth:.2f}")
            return average_depth
        except ZeroDivisionError:
            logger.debug("Length of BedRecord is 0, returning average depth of coverage as 0.0")
            return 0.0
