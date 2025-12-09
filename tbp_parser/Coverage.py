import subprocess
import globals as globals_
import pandas as pd
import os

class Coverage:
    """
    This class creates the CDPH coverage report.
    It has several functions: 
        - calculate_depth: uses samtools to calculate the breadth of coverage and average depth for a given region
        - get_coverage: creates the initial coverage report.
        - reformat_coverage: adds a warning column if a deletion was identified
            within the gene.
    """

    def __init__(self, logger, input_bam, output_prefix, coverage_regions):
        """ Initalizes the Coverage class

        Args:
            logger (logging.getlogger() object): Object that handles logging
            input_bam (File): BAM file of sample to be analyzed (aligned to H37Rv)
            output_prefix (String): Prefix for all output
            coverage_regions (File): Bed file of regions to be examined for coverage
        """

        self.logger = logger
        self.input_bam = input_bam
        self.output_prefix = output_prefix
        self.coverage_regions = coverage_regions
        
        # extract chromosome name from BAM file -- sometimes this can be different depending on reference used
        command = "samtools idxstats {} | cut -f 1 | head -1".format(self.input_bam)
        self.chromosome = subprocess.check_output(command, shell=True).decode("utf-8").strip()

    def calculate_depth(self, line) -> tuple[str, float, float]:
        """Uses samtools to calculate the breadth of coverage and average depth for a given region

        Args:
            line (str): A line from a bed file listing regions of interest

        Returns:
            tuple[str, float, float]: the name of the region, the percentage of the region over the minimum depth, and the average depth of the region
        """        
        # parse out the coordinates and gene from each line in the bed file -- 
        #  assuming 1-based indexing and formatting based on TBDB.bed from TBProfiler
        start = line[1]
        end = line[2]
        gene = line[4]

        depth_command = "samtools depth -a -J -r \"{}:{}-{}\" {} > {}-depths.txt".format(self.chromosome, start, end, self.input_bam, gene)
        subprocess.Popen(depth_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

        # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
        #  above min_depth, then wc -l counts them all
        breadth_coverage_command = "awk -F '\t' '{if ($3 >= " + str(globals_.MIN_DEPTH) + ") print;}' " + gene + "-depths.txt | wc -l"
        depth = subprocess.Popen(breadth_coverage_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
        
        # get the percentage of the region that is above the minimum depth; add one to gene length to compensate for subtraction (1-based indexing)
        breadth_of_coverage = (int(depth) / (int(end) - int(start) + 1)) * 100
        
        # add together the depth of all the sites, then divide by the number of sites to get the average depth
        average_depth_command = "awk -F '\t' '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }' " + gene + "-depths.txt"
        average_depth = subprocess.Popen(average_depth_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
        
        self.logger.debug("COV:calculate_depth:The average breadth of coverage for this gene ({}) is {}".format(gene, breadth_of_coverage))
        self.logger.debug("COV:calculate_depth:The average depth for this gene ({}) is {}".format(gene, float(average_depth)))

        # remove gene + "-depths.txt file to save space
        os.remove(gene + "-depths.txt")

        return gene, breadth_of_coverage, float(average_depth)

    def get_coverage(self) -> None:
        """ Iterates through a bedfile and adds the breadth of coverage to the global variable "COVERAGE_DICTIONARY" 
        and the average depth to "AVERAGE_LOCI_COVERAGE"

        Args: 
            None

        Returns:
            None
        """
        self.logger.debug("COV:get_coverage:The chromosome name was collected during class initialization: {}".format(self.chromosome))
        self.logger.debug("COV:get_coverage:Now calculating breadth of coverage and average coverage for each gene in the {} file".format(self.coverage_regions))

        with open(self.coverage_regions, "r") as bedfile_fh:
            for line in bedfile_fh:
                line = line.split("\t")
                gene, coverage, average_depth = self.calculate_depth(line)

                # rename genes to match CDPH nomenclature -- 
                # TO-DO: make sure that this is acceptable for everyone -- maybe make this configurable????  
                # a global rename dictionary? 
                if gene == "mmpR5":
                    gene = "Rv0678"
                elif gene == "fbiD":
                    gene = "Rv2983"

                globals_.COVERAGE_DICTIONARY[gene] = coverage
                globals_.AVERAGE_LOCI_COVERAGE[gene] = average_depth

        self.logger.info("COV:get_coverage:Coverage dictionaries of length {} and {} have been created\n".format(len(globals_.COVERAGE_DICTIONARY), len(globals_.AVERAGE_LOCI_COVERAGE)))

    def create_coverage_report(self) -> None:
        """
        This function reformats the coverage and average loci coverage dictionaries into
        a CSV file and adds a deletion warning if a QC-passing deletion was identified
        for the region in the Laboratorian report.
        """
        # TO-DO: enable customizable column headers!!!!
        DF_COVERAGE = pd.DataFrame(columns=["Gene", "Percent_Coverage", "Average_Locus_Coverage", "Warning"])

        self.logger.debug("COV:create_coverage_report:Now iterating through each gene in the breadth of coverage dictionary")
        for gene, percent_coverage in globals_.COVERAGE_DICTIONARY.items():            
            average_coverage = globals_.AVERAGE_LOCI_COVERAGE[gene]
            warning = ""

            try:
                # get a list of the nucleotide mutations for the given gene
                nucleotide_mutations = globals_.DF_LABORATORIAN["tbprofiler_variant_substitution_nt"][globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]
                for mutation in nucleotide_mutations:
                    if "del" in mutation and mutation not in globals_.MUTATION_FAIL_LIST:
                        warning = "Deletion identified"
                        if float(percent_coverage) == 100:
                            warning = "Deletion identified (upstream)"
                        self.logger.debug("COV:create_coverage_report:A deletion warning is being added to a gene ({}) with {}%% coverage: {}".format(gene, percent_coverage, warning))
            except:
                self.logger.error("An expected gene ({}) was not found in laboratorian report.\nSomething may have gone wrong.".format(gene))

            # prevent concatenation warnings if the dataframe is empty
            if len(DF_COVERAGE) == 0:
                DF_COVERAGE = pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, 
                                            "Average_Locus_Coverage": average_coverage, "Warning": warning}, 
                                           index=[0])
            else:
                DF_COVERAGE = pd.concat([DF_COVERAGE, pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, 
                                                                   "Average_Locus_Coverage": average_coverage, "Warning": warning}, 
                                                                   index=[0])], ignore_index=True)

        if globals_.TNGS: 
            # TO-DO: MAKE CONFIGURABLE
            DF_COVERAGE.rename(columns={"Percent_Coverage": "Coverage_Breadth_reportableQC_region", "Warning": "QC_Warning"}, inplace=True)

        DF_COVERAGE.to_csv(self.output_prefix + ".coverage_report.csv", index=False)
        self.logger.info("COV:create_coverage_report:The coverage report has been saved to {}\n".format(self.output_prefix + ".coverage_report.csv"))