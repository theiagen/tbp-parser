import subprocess
import globals as globals_
import pandas as pd
import os

class Coverage:
    """A class representing the coverage report and coverage statistics
    
    Attributes:
        logger (logging.getlogger() object): Object that handles logging
        input_bam (str): path to BAM file of sample to be analyzed (aligned to H37Rv)
        OUTPUT_PREFIX (str): Prefix for all output
        coverage_regions (str): path to BED file of regions to be examined for coverage
    
    Methods:
        calculate_depth(line: str, MIN_DEPTH: int) -> tuple[str, float, float]:
            uses samtools to calculate the breadth of coverage and average depth for a given region
        
        get_coverage(MIN_PERCENT_COVERAGE: float, MIN_DEPTH: int) -> tuple[dict[str, float], dict[str, float], list[str]]:
            iterates through a bedfile and adds the breadth of coverage to the coverage dictionary and the average depth to the average loci coverage dictionary
        
        create_coverage_report(COVERAGE_DICTIONARY: dict[str, float], AVERAGE_LOCI_COVERAGE: dict[str, float], GENES_WITH_VALID_DELETIONS: list[str], TNGS: bool) -> None:
            reformats the coverage and average loci coverage dictionaries into a CSV file and adds a deletion warning if a QC-passing deletion was identified for the region
    """
    def __init__(self, logger, input_bam, OUTPUT_PREFIX, coverage_regions) -> None:
        """ Initalizes the Coverage class

        Args:
            logger (logging.getlogger() object): Object that handles logging
            input_bam (str): path to BAM file of sample to be analyzed (aligned to H37Rv)
            OUTPUT_PREFIX (str): Prefix for all output
            coverage_regions (str): path to BED file of regions to be examined for coverage
        """

        self.logger = logger
        self.input_bam = input_bam
        self.OUTPUT_PREFIX = OUTPUT_PREFIX
        self.coverage_regions = coverage_regions
        
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
        
        region_length = int(end) - int(start) + 1  #1-based indexing
        breadth_of_coverage = (positions_above_min / region_length) * 100
        average_depth = (depth_sum / region_length) if region_length > 0 else 0
        
        self.logger.debug("COV:calculate_depth:The average breadth of coverage for this gene ({}) is {} and the average depth is {}".format(gene, breadth_of_coverage, average_depth))
        
        os.remove(gene + "-depths.txt")

        return gene, float(breadth_of_coverage), float(average_depth)

    def get_coverage(self, MIN_PERCENT_COVERAGE: float, MIN_DEPTH: int) -> tuple[dict[str, float], dict[str, float], list[str]]:
        """ Iterates through a bedfile and adds the breadth of coverage to the global variable "COVERAGE_DICTIONARY" 
        and the average depth to "AVERAGE_LOCI_COVERAGE"

        Args: 
            MIN_PERCENT_COVERAGE (float): The minimum percent breadth of coverage threshold for a gene/locus
            MIN_DEPTH (int): The minimum depth threshold for a gene/locus

        Returns:
            tuple[dict, dict, list]: A tuple containing the coverage dictionary, average loci coverage dictionary, and a list of genes below the MIN_PERCENT_COVERAGE
        """
        self.logger.debug("COV:get_coverage:The chromosome name was collected during class initialization: {}".format(self.chromosome))
        self.logger.debug("COV:get_coverage:Now calculating breadth of coverage and average coverage for each gene in the {} file".format(self.coverage_regions))

        COVERAGE_DICTIONARY = {}
        AVERAGE_LOCI_COVERAGE = {}

        with open(self.coverage_regions, "r") as bedfile_fh:
            for line in bedfile_fh:
                line = line.split("\t")
                gene, breadth_of_coverage, average_depth = self.calculate_depth(line, MIN_DEPTH)

                COVERAGE_DICTIONARY[gene] = breadth_of_coverage
                AVERAGE_LOCI_COVERAGE[gene] = average_depth
                
        # add to low depth of coverage list if below the breadth of coverage threshold (MIN_PERCENT_COVERAGE * 100)
        LOW_DEPTH_OF_COVERAGE_LIST = [gene for gene, coverage in COVERAGE_DICTIONARY.items() if coverage < (MIN_PERCENT_COVERAGE * 100)]

        self.logger.info("COV:get_coverage:Coverage dictionaries of length {} and {} have been created".format(len(COVERAGE_DICTIONARY), len(AVERAGE_LOCI_COVERAGE)))
        self.logger.info("COV:get_coverage:The following genes (total: {}) have coverage below the {}% breadth of coverage threshold: {}\n".format(len(LOW_DEPTH_OF_COVERAGE_LIST), (MIN_PERCENT_COVERAGE * 100), LOW_DEPTH_OF_COVERAGE_LIST))

        return COVERAGE_DICTIONARY, AVERAGE_LOCI_COVERAGE, LOW_DEPTH_OF_COVERAGE_LIST

    def create_coverage_report(self, COVERAGE_DICTIONARY: dict[str, float], AVERAGE_LOCI_COVERAGE: dict[str, float], GENES_WITH_VALID_DELETIONS: list[str], TNGS: bool) -> None:
        """
        This function reformats the coverage and average loci coverage dictionaries into
        a CSV file and adds a deletion warning if a QC-passing deletion was identified
        for the region in the Laboratorian report.
        """
        DF_COVERAGE = pd.DataFrame(columns=["Gene", "Percent_Coverage", "Average_Locus_Coverage", "Warning"])

        self.logger.debug("COV:create_coverage_report:Now iterating through each gene in the breadth of coverage dictionary")
        for gene, percent_coverage in COVERAGE_DICTIONARY.items():            
            average_coverage = AVERAGE_LOCI_COVERAGE[gene]
            warning = ""

            if gene in GENES_WITH_VALID_DELETIONS:
                warning = "Deletion identified"
                if percent_coverage == 100: 
                    # if the coverage is at 100% but a deletion was identified, it's likely upstream
                    warning = "Deletion identified (upstream)"

            # prevent concatenation warnings if the dataframe is empty
            if len(DF_COVERAGE) == 0:
                DF_COVERAGE = pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, 
                                            "Average_Locus_Coverage": average_coverage, "Warning": warning}, 
                                           index=[0])
            else:
                DF_COVERAGE = pd.concat([DF_COVERAGE, pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, 
                                                                   "Average_Locus_Coverage": average_coverage, "Warning": warning}, 
                                                                   index=[0])], ignore_index=True)

        if TNGS: 
            DF_COVERAGE.rename(columns={"Percent_Coverage": "Coverage_Breadth_reportableQC_region", "Warning": "QC_Warning"}, inplace=True)

        DF_COVERAGE.to_csv(self.OUTPUT_PREFIX + ".coverage_report.csv", float_format="%.2f", index=False)