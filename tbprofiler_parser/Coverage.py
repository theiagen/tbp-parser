import subprocess
import globals
import pandas as pd
import importlib_resources

class Coverage:
  def __init__(self, logger, input_bam, output_prefix):
    self.logger = logger
    self.input_bam = input_bam
    self.output_prefix = output_prefix
  
  def calculate_coverage(self):
    """
    This function calculates the coverage for each gene in the tbdb.bed file.
    """
    
    self.logger.info("Within calculate_coverage function")
    self.logger.info("Collecting the chromosome name from the BAM file")
    
    command = "samtools idxstats {} | cut -f 1 | head -1".format(self.input_bam)
    CHROMOSOME = subprocess.check_output(command, shell=True).decode("utf-8").strip()
    #self.logger.debug("The BAM file Chromosome: {}".format(CHROMOSOME))    
    
    with open(importlib_resources.files(__name__) / "../data/tbdb.bed", "r") as bedfile_fh:
      for line in bedfile_fh:
        line = line.split("\t")
        start = line[1]
        end = line[2]
        gene = line[4]
        
        # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
        #  above min_depth, then wc -l counts them all
        command = "samtools depth -J -r \"" + CHROMOSOME + ":" + start + "-" + end + "\" " + self.input_bam + " | awk -F '\t' '{if ($3 >= " + str(globals.MIN_DEPTH) + ") print;}' | wc -l"
        depth = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
        
        # get coverage for region in bed file based on depth
        #  add one to gene length to compensate for subtraction
        coverage = (int(depth) / (int(end) - int(start) + 1)) * 100
        #self.logger.debug("Gene: {} Coverage: {}".format(gene, coverage))
        
        globals.COVERAGE_DICTIONARY[gene] = coverage
    
    # rename some genes to match CDPH preference
    globals.COVERAGE_DICTIONARY["Rv0678"] = globals.COVERAGE_DICTIONARY["mmpR5"]
    globals.COVERAGE_DICTIONARY["Rv2983"] = globals.COVERAGE_DICTIONARY["fbiD"]
    
    #self.logger.debug("Coverage dictionary: ")
    #self.logger.debug(globals.COVERAGE_DICTIONARY)
    
  def reformat_coverage(self):
    """
    This function reformats the coverage dictionary into a CSV file
    and includes the warning field.
    """
    self.logger.info("Within reformat_coverage function")
    
    DF_COVERAGE = pd.DataFrame(columns=["Gene", "Percent_Coverage", "Warning"])

    for gene, percent_coverage in globals.COVERAGE_DICTIONARY.items():
      warning = ""
      
      try:
        for mutation_type_nucleotide in globals.DF_LABORATORIAN["tbprofiler_variant_substitution_nt"][globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]:
          if "del" in mutation_type_nucleotide:
            warning = "Deletion identified"
            if float(percent_coverage) == 100:
              warning = "Deletion identified (upstream)"
      except:
        self.logger.debug("Gene {} not found in lab report".format(gene))
      
      
      self.logger.debug("Gene: {} Percent Coverage: {} Warning: {}".format(gene, percent_coverage, warning))
      DF_COVERAGE = pd.concat([DF_COVERAGE, pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, "Warning": warning}, index=[0])], ignore_index=True)

    DF_COVERAGE.to_csv(self.output_prefix + ".percent_gene_coverage.csv", index=False)