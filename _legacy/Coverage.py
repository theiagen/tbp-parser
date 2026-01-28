import subprocess
import pandas as pd
import globals as globals_
import copy
import os

class Coverage:
    """A class representing the coverage report and coverage statistics
    
    Attributes:
        logger (logging.getlogger() object): Object that handles logging
        input_bam (str): path to BAM file of sample to be analyzed (aligned to H37Rv)
        OUTPUT_PREFIX (str): Prefix for all output
        tbdb_bed (str): path to BED file of regions to be examined for coverage
        TNGS_REGIONS (dict[str, list[int] | dict[str, list[int]]]): A dictionary containing the specific tNGS regions for each gene (to account for non-overlapping primers in the same gene).
        USE_ERR_AS_BRR (bool): Whether to use essential region coverage as breadth of coverage for QC
        err_bed (str): path to BED file of essential regions to be examined for coverage
        
        err_regions (dict[str, list[int]]): A dictionary containing the start and end positions of the regions in the err_bed file
        err_coverage_dictionary (dict[str, float]): A dictionary containing the breadth of coverage for the regions in the err_bed file
        brr_coverage_dictionary (dict[str, float]): A dictionary containing the breadth of coverage for the regions in the tbdb_bed file
        average_loci_coverage (dict[str, float]): A dictionary containing the average depth of coverage for the regions in the tbdb_bed file
        chromosome (str): The chromosome name extracted from the BAM file
    
    Methods:
        calculate_depth(line: str, MIN_DEPTH: int) -> tuple[str, float, float]:
            uses samtools to calculate the breadth of coverage and average depth for a given region
        
        get_coverage(MIN_PERCENT_COVERAGE: float, MIN_DEPTH: int) -> tuple[dict[str, float], list[str]]:
            iterates through a bedfile and adds the breadth of coverage to the coverage dictionary and the average depth to the average loci coverage dictionary
        
        create_coverage_report(COVERAGE_DICTIONARY: dict[str, float], AVERAGE_LOCI_COVERAGE: dict[str, float], GENES_WITH_VALID_DELETIONS: dict[str, list[int]], TNGS: bool) -> None:
            reformats the coverage and average loci coverage dictionaries into a CSV file and adds a deletion warning if a QC-passing deletion was identified for the region
    """
    
    def __init__(self, logger, input_bam: str, OUTPUT_PREFIX: str, tbdb_bed: str, TNGS_REGIONS: dict[str, list[int] | dict[str, list[int]]], USE_ERR_AS_BRR: bool, err_bed: str = None) -> None:
        """ Initalizes the Coverage class

        Args:
            logger (logging.getlogger() object): Object that handles logging
            input_bam (str): path to BAM file of sample to be analyzed (aligned to H37Rv)
            OUTPUT_PREFIX (str): Prefix for all output
            tbdb_bed (str): path to BED file of regions to be examined for coverage\
            TNGS_REGIONS (dict[str, list[int] | dict[str, list[int]]]): A dictionary containing the specific tNGS regions for each gene (to account for non-overlapping primers in the same gene).
        """

        self.logger = logger
        self.input_bam = input_bam
        self.OUTPUT_PREFIX = OUTPUT_PREFIX
        self.tbdb_bed = tbdb_bed
        
        self.TNGS_REGIONS = TNGS_REGIONS
        self.USE_ERR_AS_BRR = USE_ERR_AS_BRR
        self.err_bed = err_bed
        
        self.err_regions = {}
        self.err_coverage_dictionary = {}
        self.brr_coverage_dictionary = {}
        self.average_loci_coverage = {}
        
        # extract chromosome name from BAM file -- sometimes this can be different depending on reference used
        command = "samtools idxstats {} | cut -f 1 | head -1".format(self.input_bam)
        self.chromosome = subprocess.check_output(command, shell=True).decode("utf-8").strip()

    def calculate_depth(self, line: str, MIN_DEPTH: int) -> tuple[str, float, float]:
        """Uses samtools to calculate the breadth of coverage and average depth for a given region

        Args:
            line (str): A line from a bed file listing regions of interest
            MIN_DEPTH (int): The minimum depth threshold for a gene/locus

        Returns:
            tuple[str, float, float]: the name of the region, the percentage of the region over the minimum depth, and the average depth of the region
        """        
        # parse out the coordinates and gene from each line in the bed file -- 
        #  assuming 1-based indexing and formatting based on TBDB.bed from TBProfiler
        start = line[1]
        end = line[2]
        gene = line[4]

        depth_command = "samtools depth -a -J -r \"{}:{}-{}\" {} > {}-depths.txt".format(self.chromosome, start, end, self.input_bam, gene)
        subprocess.run(depth_command, shell=True, check=True)

        # use python to calculate the statistics instead of using awk to prevent subprocess confusion
        total_positions = 0
        positions_above_min = 0
        depth_sum = 0

        with open(gene + "-depths.txt", "r") as depths_fh:
            for line in depths_fh:
                cols = line.strip().split()
                depth = int(cols[2])
                total_positions += 1
                depth_sum += depth
                if depth >= MIN_DEPTH:
                    positions_above_min += 1
        
        region_length = int(end) - int(start) + 1  # 1-based indexing
        
        if region_length > 0:
            breadth_of_coverage = (positions_above_min / region_length) * 100
            average_depth = depth_sum / region_length
        else:
            breadth_of_coverage = 0
            average_depth = 0
            
        self.logger.debug("COV:calculate_depth:The average breadth of coverage for this gene ({}) is {} and the average depth is {}".format(gene, breadth_of_coverage, average_depth))
        
        os.remove(gene + "-depths.txt")

        return gene, float(breadth_of_coverage), float(average_depth)

    def get_coverage(self, MIN_PERCENT_COVERAGE: float, MIN_DEPTH: int, TNGS: bool) -> tuple[dict[str, float], list[str]]:
        """ Iterates through a bedfile and adds the breadth of coverage to the coverage_dictionary
            and the average depth to the average_loci_coverage dictionary. If an err_bed file was provided,
            it also calculates coverage for those regions and adds them to the err_coverage object dictionary.

        Args: 
            MIN_PERCENT_COVERAGE (float): The minimum percent breadth of coverage threshold for a gene/locus
            MIN_DEPTH (int): The minimum depth threshold for a gene/locus

        Returns:
            tuple[dict, list]: A tuple containing the coverage dictionary and a list of genes below the MIN_PERCENT_COVERAGE
        """
        self.logger.debug("COV:get_coverage:The chromosome name was collected during class initialization: {}".format(self.chromosome))
        self.logger.debug("COV:get_coverage:Now calculating breadth of coverage and average coverage for each gene in the {} file".format(self.tbdb_bed))
        
        low_depth_of_coverage_list = []

        with open(self.tbdb_bed, "r") as bedfile_fh:
            for line in bedfile_fh:
                line = line.strip().split("\t")
                gene, breadth_of_coverage, average_depth = self.calculate_depth(line, MIN_DEPTH)

                self.brr_coverage_dictionary[gene] = breadth_of_coverage
                self.average_loci_coverage[gene] = average_depth

        if self.err_bed is not None:
            with open(self.err_bed, "r") as err_bed_fh:
                for line in err_bed_fh:
                    line = line.strip().split("\t")
                    start = line[1]
                    end = line[2]
                    gene, breadth_of_coverage, average_depth = self.calculate_depth(line, MIN_DEPTH)

                    self.err_regions[gene] = [int(start), int(end)]  # store the start and end positions of the essential region
                    self.err_coverage_dictionary[gene] = breadth_of_coverage
                    
        if self.USE_ERR_AS_BRR and self.err_bed is not None:
            coverage_dictionary = copy.deepcopy(self.err_coverage_dictionary)
        else:
            coverage_dictionary = copy.deepcopy(self.brr_coverage_dictionary)
    
        # add to low depth of coverage list if below the breadth of coverage threshold (MIN_PERCENT_COVERAGE * 100)
        # combine split primers so if one fails, gene fails rpoB_1 = 0 but rpob_2 = 100 -> rpoB is added to this list
        if TNGS:
            for parent_gene, children in self.TNGS_REGIONS.items():
                if isinstance(children, dict):
                    # 7.1.3.1 - breadth of coverage QC fails if at least one segment does not meet QC thresholds (use the minimum here)
                    combined_coverage = min(coverage_dictionary[region] for region in children)
                    if combined_coverage < (MIN_PERCENT_COVERAGE * 100):
                        low_depth_of_coverage_list.append(parent_gene)
                        
        # this will contain split primers if tNGS, but they will be reported as part of the parent gene above
        low_depth_of_coverage_list.extend([gene for gene, coverage in coverage_dictionary.items() if coverage < (MIN_PERCENT_COVERAGE * 100)])
    
        self.logger.info("COV:get_coverage:Coverage dictionaries of length {} and {} have been created".format(len(coverage_dictionary), len(self.average_loci_coverage)))
        self.logger.info("COV:get_coverage:The following genes (total: {}) have coverage below the {}% breadth of coverage threshold: {}\n".format(len(low_depth_of_coverage_list), (MIN_PERCENT_COVERAGE * 100), low_depth_of_coverage_list))

        return coverage_dictionary, low_depth_of_coverage_list

    def create_coverage_report(self, SAMPLE_NAME: str, GENES_WITH_VALID_DELETIONS: dict[str, list[int]]) -> None:
        """
        This function reformats the coverage and average loci coverage dictionaries into
        a CSV file and adds a deletion warning if a QC-passing deletion was identified
        for the region in the Laboratorian report.
        """
        self.logger.debug("COV:create_coverage_report:Now iterating through each gene in the breadth of coverage dictionary")
        
        coverage_report = []
        
        for gene, percent_coverage in self.brr_coverage_dictionary.items():            
            average_coverage = self.average_loci_coverage[gene]
            warning = ""

            if gene in GENES_WITH_VALID_DELETIONS.keys():
                if self.USE_ERR_AS_BRR and gene in self.err_regions.keys():
                    for positions in  GENES_WITH_VALID_DELETIONS[gene]:
                        if globals_.is_mutation_within_range(list(positions), self.err_regions[gene]):
                            warning = "Deletion identified"
                            if percent_coverage == 100: 
                                warning = "Deletion identified (upstream)"
                            self.logger.info("COV:create_coverage_report:Deletion positions {} for gene {} are within the essential region {}; adding a warning".format(positions, gene, self.err_regions[gene]))
                            # at least one deletion position is within the essential region; no need to check further
                            break
                        else:
                            self.logger.info("COV:create_coverage_report:Deletion positions {} for gene {} are outside the essential region {}; not adding a warning".format(positions, gene, self.err_regions[gene]))
                else:                    
                    warning = "Deletion identified"
                    if percent_coverage == 100: 
                        # if the coverage is at 100% and a deletion was identified, it's likely upstream
                        warning = "Deletion identified (upstream)"

            if len(self.err_coverage_dictionary) == 0:
                coverage_report.append({
                    "Sample": SAMPLE_NAME,
                    "Gene": gene,
                    "Percent_Coverage": percent_coverage,
                    "Average_Locus_Coverage": average_coverage,
                    "Warning": warning
                })
            else:
                err_coverage = self.err_coverage_dictionary.get(gene, "N/A")
                
                coverage_report.append({
                    "Sample": SAMPLE_NAME,
                    "Gene": gene,
                    "Percent_Coverage": percent_coverage,
                    "Average_Locus_Coverage": average_coverage,
                    "Essential_Regions_Percent_Coverage": err_coverage,
                    "Warning": warning
                })

        DF_COVERAGE = pd.DataFrame(coverage_report)
        DF_COVERAGE.to_csv(self.OUTPUT_PREFIX + ".coverage_report.csv", float_format="%.2f", index=False)
