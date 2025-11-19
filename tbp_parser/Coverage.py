import subprocess
import globals as globals_
import pandas as pd

class Coverage:
  """
  This class creates the CDPH coverage report.
  It has two functions: 
    - calculate_coverage: creates the initial coverage report.
    - reformat_coverage: adds a warning column if a deletion was identified
      within the gene.
  """
  
  def __init__(self, logger, input_bam, output_prefix, coverage_regions, tngs_expert_regions):
    """ Initalizes the Coverage class

    Args:
      logger (logging.getlogger() object): Object that handles logging
      input_bam (File): BAM file of sample to be analyzed (aligned to H37Rv)
      output_prefix (String): Prefix for all output
      coverage_regions (File): Bed file of regions to be examined for coverage
      tngs_expert_regions (File): Bed file of the expert rule regions to be examined for coverage in tNGS data
    """
    
    self.logger = logger
    self.input_bam = input_bam
    self.output_prefix = output_prefix
    self.coverage_regions = coverage_regions
    self.tngs_expert_regions = tngs_expert_regions
    self.tngs_expert_regions_coverage = {}
    
    command = "samtools idxstats {} | cut -f 1 | head -1".format(self.input_bam)
    self.chromosome = subprocess.check_output(command, shell=True).decode("utf-8").strip()
  
  def calculate_depth(self, line):
    """ Uses samtools to calculate average breadth of coverage

    Args:
      line (String): A line from a bed file listing regions of interest

    Returns:
      String gene: name of the region
      Float coverage: average coverage of the region over minimum depth
    """
    
    # parse out the coordinates and gene from each line in the bed file
    start = line[1]
    end = line[2]
    gene = line[4]
    
    # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
    #  above min_depth, then wc -l counts them all
    command = "samtools depth -r \"" + self.chromosome + ":" + start + "-" + end + "\" " + self.input_bam + " | awk -F '\t' '{if ($3 >= " + str(globals_.MIN_DEPTH) + ") print;}' | wc -l"
    self.logger.debug("COV:Now running " + command)
    depth = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
    
    # get coverage for region in bed file based on depth
    #  add one to gene length to compensate for subtraction
    coverage = (int(depth) / (int(end) - int(start) + 1)) * 100
    self.logger.debug("COV:The coverage for this gene ({}) is {}".format(gene, coverage))
    return gene, coverage
  
  def calculate_average_depth(self, line):
    """ Uses samtools to calculate average coverage of a locus

    Args:
      line (String): A line from a bed file listing regions of interest

    Returns:
      String gene: name of the region
      Float coverage: average coverage of the region over minimum depth
    """
    
    # parse out the coordinates and gene from each line in the bed file
    start = line[1]
    end = line[2]
    gene = line[4]
    
    # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
    #  above min_depth, then wc -l counts them all
    command = "samtools depth -a -J -r \"" + self.chromosome + ":" + start + "-" + end + "\" " + self.input_bam + " | awk -F '\t' '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }'"
    self.logger.debug("COV:Now running " + command)
    average_depth = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
    
    self.logger.debug("COV:The average_depth for this locus ({}) is {}".format(gene, average_depth))
    return gene, average_depth
  
  def calculate_coverage(self):
    """ Iterates through a bedfile and adds average breadth of coverage to global variable "COVERAGE_DICTIONARY"
    
    Args: 
      None
        
    Returns:
      None
    """    

    self.logger.info("COV:Within the Coverage class calculate_coverage function")
    self.logger.debug("COV:The chromosome name was collected during class initialization: {}".format(self.chromosome))
    
    self.logger.debug("COV:The BAM file chromosome name is: {}".format(self.chromosome))    
    
    with open(self.coverage_regions, "r") as bedfile_fh:
      self.logger.debug("COV:Now calculating coverage for each gene in the {} file".format(self.coverage_regions))
      for line in bedfile_fh:
        line = line.split("\t")
        gene, coverage = self.calculate_depth(line)
        
        globals_.COVERAGE_DICTIONARY[gene] = coverage
    
    # rename some genes to match CDPH nomenclature
    if "mmpR5" in globals_.COVERAGE_DICTIONARY.keys():
      globals_.COVERAGE_DICTIONARY["Rv0678"] = globals_.COVERAGE_DICTIONARY["mmpR5"]
    if "fbiD" in globals_.COVERAGE_DICTIONARY.keys():
      globals_.COVERAGE_DICTIONARY["Rv2983"] = globals_.COVERAGE_DICTIONARY["fbiD"]
    
    self.logger.info("COV:Initial coverage report created, now exiting function\n")
   
  def calculate_r_expert_rule_regions_coverage(self):
    
    """
    This function calculates the breadth of coverage over the ranges that encompass
    R mutations and expert rule regions; intended for supervisory review in the case
    of low breadth of coverage and should not be used as a QC threshold.
    """
    self.logger.info("COV:Within the Coverage class calculate_r_expert_rule_regions_coverage function")
    self.logger.debug("COV:Now calculating coverage with the tngs-expert-rule-regions.bed file")
    
    with open(self.tngs_expert_regions, "r") as bedfile_fh:
      self.logger.debug("COV:Now calculating coverage for each gene in the {} file".format(self.tngs_expert_regions))
      for line in bedfile_fh:
        line = line.split("\t")
        gene, coverage = self.calculate_depth(line)
        
        self.tngs_expert_regions_coverage[gene] = coverage
          
    # rename some genes to match CDPH nomenclature
    if "mmpR5" in globals_.COVERAGE_DICTIONARY.keys():
       self.tngs_expert_regions_coverage["Rv0678"] =  self.tngs_expert_regions_coverage["mmpR5"]
        
    self.logger.info("COV:Expert regions coverage dictionary created, now exiting function\n")
    
  def reformat_coverage(self):
    """
    This function reformats the coverage dictionary into a CSV file and 
    includes a deletion warning field as determined by the Laboratorian report.
    """
    self.logger.info("COV:Within the Coverage class reformat_coverage function")
    DF_COVERAGE = pd.DataFrame(columns=["Gene", "Percent_Coverage", "Warning"])

    self.logger.debug("COV:Now iterating through each gene in the inital coverage report")
    for gene, percent_coverage in globals_.COVERAGE_DICTIONARY.items():
      warning = ""
      
      try:
        for mutation_type_nucleotide in globals_.DF_LABORATORIAN["tbprofiler_variant_substitution_nt"][globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]:
          if "del" in mutation_type_nucleotide and mutation_type_nucleotide not in globals_.MUTATION_FAIL_LIST:
            warning = "Deletion identified"
            if float(percent_coverage) == 100:
              warning = "Deletion identified (upstream)"
            self.logger.debug("COV:A deletion warning is being added to a gene ({}) with {}%% coverage: {}".format(gene, percent_coverage, warning))
      except:
        self.logger.error("An expected gene ({}) was not found in laboratorian report.\nSomething may have gone wrong.".format(gene))
      
      if len(DF_COVERAGE) == 0:
        DF_COVERAGE = pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, "Warning": warning}, index=[0])
      else:
        DF_COVERAGE = pd.concat([DF_COVERAGE, pd.DataFrame({"Gene": gene, "Percent_Coverage": percent_coverage, "Warning": warning}, index=[0])], ignore_index=True)

    if globals_.TNGS:
      self.logger.debug("COV:Merging the tNGS expert rule regions coverage with the initial coverage report and renaming columns")
      df_tngs_expert_regions_coverage = pd.DataFrame(self.tngs_expert_regions_coverage, index=[0]).T.reset_index().rename(columns={"index": "Gene", 0: "Coverage_Breadth_R_expert-rule_region"})
      DF_COVERAGE = pd.merge(DF_COVERAGE, df_tngs_expert_regions_coverage, on="Gene", how="outer")
      DF_COVERAGE.rename(columns={"Percent_Coverage": "Coverage_Breadth_reportableQC_region", "Warning": "QC_Warning"}, inplace=True)
      
      # add average_locus_coverage
      with open(self.coverage_regions, "r") as bedfile_fh:
        self.logger.debug("COV:Now calculating locus coverage for each gene in the {} file".format(self.coverage_regions))
        average_depths = {}
        for line in bedfile_fh:
          line = line.split("\t")
          gene, average_depth = self.calculate_depth(line)
          average_depths[gene] = average_depth
          df_average_depths = pd.DataFrame(average_depths, index=[0]).T.reset_index().rename(columns={"index": "Gene", 0: "Average_Locus_Coverage"})
          DF_COVERAGE = pd.merge(DF_COVERAGE, df_average_depths, on="Gene", how="outer")      

    DF_COVERAGE.to_csv(self.output_prefix + ".percent_gene_coverage.csv", index=False)
    self.logger.info("COV:Coverage report reformatted and saved to {}\n".format(self.output_prefix + ".percent_gene_coverage.csv"))