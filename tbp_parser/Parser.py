import logging
import subprocess
import globals as globals_
import sys
import yaml
from Coverage import Coverage
from Laboratorian import Laboratorian
from Looker import Looker
from LIMS import LIMS


class Parser:
  """This class orchestrates the different modules within the tbp_parser tool."""
  
  def __init__(self, options):
    logging.basicConfig(encoding='utf-8', level=logging.ERROR, stream=sys.stderr)
    self.logger = logging.getLogger(__name__)
    self.input_json = options.input_json
    self.input_bam = options.input_bam
    self.config = options.config
    self.tngs = options.tngs
    self.verbose = options.verbose
    self.debug = options.debug
    self.output_prefix = options.output_prefix
    self.coverage_regions = options.coverage_regions
    self.tngs_expert_regions = options.tngs_expert_regions
    self.add_cs_lims = options.add_cs_lims
    globals_.MIN_DEPTH = options.min_depth
    globals_.COVERAGE_THRESHOLD = options.min_percent_coverage
    globals_.SEQUENCING_METHOD = options.sequencing_method
    globals_.MIN_READ_SUPPORT = options.min_read_support
    globals_.MIN_FREQUENCY = options.min_frequency
    globals_.RRS_FREQUENCY = options.rrs_frequency
    globals_.RRS_READ_SUPPORT = options.rrs_read_support
    globals_.RRL_FREQUENCY = options.rrl_frequency
    globals_.RRL_READ_SUPPORT = options.rrl_read_support
    globals_.RPOB449_FREQUENCY = options.rpob449_frequency
    globals_.ETHA237_FREQUENCY = options.etha237_frequency
    globals_.OPERATOR = options.operator

    if self.verbose:
      self.logger.setLevel(logging.INFO)
      self.logger.info("PARSER:Verbose mode enabled")
    else:
      self.logger.setLevel(logging.ERROR)
      
    if self.debug:
      self.logger.setLevel(logging.DEBUG)
      self.logger.debug("PARSER:Debug mode enabled")

    if self.tngs:
      self.logger.info("PARSER:tNGS flag detected; adjusting outputs to reflect this")
      if (self.coverage_regions == "../data/tbdb-modified-regions.bed"):
        self.logger.debug("PARSER:Changing default coverage regions to ../data/tngs-reportable-regions.bed")
        self.coverage_regions = "../data/tngs-reportable-regions.bed"
      
      self.logger.debug("PARSER:Altering the ANTIMICROBIAL_CODE_TO_GENES dictionary to include only tNGS entries")
      globals_.ANTIMICROBIAL_CODE_TO_GENES = globals_.ANTIMICROBIAL_CODE_TO_GENES_tNGS
      
      self.logger.debug("PARSER:Altering the GENES_FOR_LIMS list to include only tNGS genes")
      globals_.GENES_FOR_LIMS = globals_.GENES_FOR_LIMS_tNGS
      
      self.logger.debug("PARSER:Setting the tNGS regions dictionary")
      globals_.TNGS_REGIONS = globals_.TNGS_REGIONS_ACTIVATED
      
    else:
      self.logger.debug("PARSER:Setting the ANTIMICROBIAL_CODE_TO_GENES dictionary to include all WGS entries")
      globals_.ANTIMICROBIAL_CODE_TO_GENES = globals_.ANTIMICROBIAL_CODE_TO_GENES_WGS
    
      self.logger.debug("PARSER:Setting the GENES_FOR_LIMS list to include all WGS genes")
      globals_.GENES_FOR_LIMS = globals_.GENES_FOR_LIMS_WGS
    
    if self.add_cs_lims:
      self.logger.info("PARSER:Adding cycloserine (CS) fields to the LIMS report")
      globals_.GENES_FOR_LIMS.extend(globals_.GENES_FOR_LIMS_CS)
      globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME.update(globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME_CS) 
      globals_.ANTIMICROBIAL_CODE_TO_GENES.update(globals_.ANTIMICROBIAL_CODE_TO_GENES_CS)
      
    if self.config != "":
      self.logger.info("PARSER:Overwriting variables with the provided config file")
      self.overwrite_variables()

  def overwrite_variables(self):
    """This function overwrites the input variables provided at runtime with those from the config file"""
    with open(self.config, "r") as config:
      settings = yaml.safe_load(config)
      
      for key, value in settings.items():
        if key.replace("self.", "") in vars(self):
          setattr(self, key.replace("self.", ""), value)
          self.logger.info("PARSER:self.{} has been overwritten with a config-specified value".format(key))
        if key.replace("globals.", "") in dir(globals_):
          setattr(globals_, key.replace("globals.", ""), value) 
          self.logger.info("PARSER:globals.{} has been overwritten with a config-specified value".format(key))
  
  def check_dependency_exists(self):
    """This function confirms that samtools is installed and available"""
    result = subprocess.run(
      ["samtools", "--version"],
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True,
    )
    
    if result.returncode != 0:
      self.logger.critical("PARSER:Error: samtools not found. Please install samtools and try again.")
      sys.exit(1)
    self.logger.info("PARSER:samtools was found! Proceeding with parsing")
  
  def run(self):
    """
    This function runs the parsing module for the tb_parser tool.
    """    
    self.logger.info("PARSER:Checking for samtools...")
    self.check_dependency_exists
    
    self.logger.info("PARSER:Creating initial coverage report")
    coverage = Coverage(self.logger, self.input_bam, self.output_prefix, self.coverage_regions, self.tngs, self.tngs_expert_regions)
    coverage.calculate_coverage()
    
    if self.tngs and self.tngs_expert_regions != "":
      self.logger.info("PARSER:Calculating the coverage for the expert rule regions")
      coverage.calculate_r_expert_rule_regions_coverage()
  
    self.logger.info("PARSER:Creating laboratorian report")
    laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix, self.tngs)
    laboratorian.create_laboratorian_report()
    
    self.logger.info("PARSER:Creating LIMS report")
    lims = LIMS(self.logger, self.input_json, self.output_prefix, self.tngs)
    lims.create_lims_report()
    
    self.logger.info("PARSER:Creating Looker report")
    looker = Looker(self.logger, self.output_prefix)
    looker.create_looker_report()
    
    self.logger.info("PARSER:Finalizing coverage report")
    coverage.reformat_coverage()
    
    self.logger.info("PARSER:Parsing completed")