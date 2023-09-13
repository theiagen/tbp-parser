import subprocess
import globals
import pandas as pd
import importlib_resources

class Coverage:
  """
  This class creates the CDPH coverage report.
  It has two functions: 
    - calculate_coverage: creates the initial coverage report.
    - reformat_coverage: adds a warning column if a deletion was identified
      within the gene.
  """
  def __init__(self, logger, input_bam, output_prefix):
    self.logger = logger
    self.input_bam = input_bam
    self.output_prefix = output_prefix
  
  def calculate_coverage(self):
    """
    This function calculates the coverage for each gene in the tbdb.bed 
    file from Jody Phelan's TBProfiler repository using `samtools depth` 
    and adds it to the global variable called "COVERAGE_DICTIONARY"
    """
    
    self.logger.info("Within the Coverage class calculate_coverage function")
    self.logger.debug("Now collecting the chromosome name from the BAM file")
    
    command = "samtools idxstats {} | cut -f 1 | head -1".format(self.input_bam)
    CHROMOSOME = subprocess.check_output(command, shell=True).decode("utf-8").strip()
    self.logger.debug("The BAM file chromosome name is: {}".format(CHROMOSOME))    
    
    with open(importlib_resources.files(__name__) / "../data/tbdb.bed", "r") as bedfile_fh:
      self.logger.debug("Now calculating coverage for each gene in the tbdb.bed file")
      for line in bedfile_fh:
        line = line.split("\t")
        start = line[1]
        end = line[2]
        gene = line[4]
        
        # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
        #  above min_depth, then wc -l counts them all
        command = "samtools depth -r \"" + CHROMOSOME + ":" + start + "-" + end + "\" " + self.input_bam + " | awk -F '\t' '{if ($3 >= " + str(globals.MIN_DEPTH) + ") print;}' | wc -l"
        self.logger.debug("Now running " + command)
        depth = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
        
        # get coverage for region in bed file based on depth
        #  add one to gene length to compensate for subtraction
        coverage = (int(depth) / (int(end) - int(start) + 1)) * 100
        self.logger.debug("The coverage for this gene ({}) is {}".format(gene, coverage))
        
        globals.COVERAGE_DICTIONARY[gene] = coverage
    
    # rename some genes to match CDPH nomenclature
    globals.COVERAGE_DICTIONARY["Rv0678"] = globals.COVERAGE_DICTIONARY["mmpR5"]
    globals.COVERAGE_DICTIONARY["Rv2983"] = globals.COVERAGE_DICTIONARY["fbiD"]
    self.logger.info("Initial coverage report created, now exiting function")
    
  def reformat_coverage(self):
    """
    This function reformats the coverage dictionary into a CSV file and 
    includes a deletion warning field as determined by the Laboratorian report.
    """
    self.logger.info("Within the Coverage class reformat_coverage function")
    DF_COVERAGE = pd.DataFrame(columns=["Gene", "Percent_Coverage", "Warning"])

    self.logger.debug("Now iterating through each gene in the inital coverage report")
    for gene, percent_coverage in globals.COVERAGE_DICTIONARY.items():
      warning = ""
      
      try:
        for mutation_type_nucleotide in globals.DF_LABORATORIAN["tbprofiler_variant_substitution_nt"][globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]:
          if "del" in mutation_type_nucleotide:
            warning = "Deletion identified"
            if float(percent_coverage) == 100:
              warning = "Deletion identified (upstream)"
            self.logger.debug("A deletion warning is being added to a gene ({}) with {}%% coverage: {}".format(gene, percent_coverage, warning))
      except:
        self.logger.error("An expected gene ({}) was not found in laboratorian report.\nSomething may have gone wrong.".format(gene))
      
      DF_COVERAGE = pd.concat([DF_COVERAGE, pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, "Warning": warning}, index=[0])], ignore_index=True)

    DF_COVERAGE.to_csv(self.output_prefix + ".percent_gene_coverage.csv", index=False)
    self.logger.info("Coverage report reformatted and saved to {}".format(self.output_prefix + ".percent_gene_coverage.csv"))