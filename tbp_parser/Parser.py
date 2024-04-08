import logging
import subprocess
import globals
import sys
from Coverage import Coverage
from Laboratorian import Laboratorian
from Looker import Looker
from LIMS import LIMS


class Parser:
  """
  This class runs the parsing module for the tb_parser tool.
  """
  def __init__(self, options):
    logging.basicConfig(encoding='utf-8', level=logging.ERROR, stream=sys.stderr)
    self.logger = logging.getLogger(__name__)
    self.input_json = options.input_json
    self.input_bam = options.input_bam
    self.tngs = options.tngs
    self.verbose = options.verbose
    self.debug = options.debug
    self.output_prefix = options.output_prefix
    self.coverage_regions = options.coverage_regions
    self.tngs_expert_regions = options.tngs_expert_regions
    globals.MIN_DEPTH = options.min_depth
    globals.COVERAGE_THRESHOLD = options.min_percent_coverage
    globals.SEQUENCING_METHOD = options.sequencing_method
    globals.MIN_READ_SUPPORT = options.min_read_support
    globals.MIN_FREQUENCY = options.min_frequency
    globals.RRS_FREQUENCY = options.rrs_frequency
    globals.RRS_READ_SUPPORT = options.rrs_read_support
    globals.RRL_FREQUENCY = options.rrl_frequency
    globals.RRL_READ_SUPPORT = options.rrl_read_support
    globals.RPOB449_FREQUENCY = options.rpob449_frequency
    globals.ETHA237_FREQUENCY = options.etha237_frequency
    globals.OPERATOR = options.operator
    
    if self.tngs:
      self.logger.info("PARSER: Deeplex + CDPH modified protocol flag detected; adjusting outputs to reflect this")
      if (self.coverage_regions == "../data/tbdb-modified-regions.bed"):
        self.logger.debug("PARSER: Changing default coverage regions to ../data/tngs-reportable-regions.bed")
        self.coverage_regions = "../data/tngs-reportable-regions.bed"
      
      self.logger.debug("PARSER: Altering the ANTIMICROBIAL_CODE_TO_GENES dictionary to include only tNGS entries")
      globals.ANTIMICROBIAL_CODE_TO_GENES = globals.ANTIMICROBIAL_CODE_TO_GENES_tNGS
      
      self.logger.debug("PARSER: Altering the GENES_FOR_LIMS list to include only tNGS genes")
      globals.GENES_FOR_LIMS = globals.GENES_FOR_LIMS_tNGS
      
      self.logger.debug("PARSER: Setting the tNGS regions dictionary")
      globals.TNGS_REGIONS = globals.TNGS_REGIONS_ACTIVATED
      
    else:
      self.logger.debug("PARSER: Setting the ANTIMICROBIAL_CODE_TO_GENES dictionary to include all WGS entries")
      globals.ANTIMICROBIAL_CODE_TO_GENES = globals.ANTIMICROBIAL_CODE_TO_GENES_WGS
    
      self.logger.debug("PARSER: Setting the GENES_FOR_LIMS list to include all WGS genes")
      globals.GENES_FOR_LIMS = globals.GENES_FOR_LIMS_WGS
    
    if self.verbose:
      self.logger.setLevel(logging.INFO)
      self.logger.info("PARSER: Verbose mode enabled")
    else:
      self.logger.setLevel(logging.ERROR)
      
    if self.debug:
      self.logger.setLevel(logging.DEBUG)
      self.logger.debug("PARSER: Debug mode enabled")
  
  def check_dependency_exists(self):
    result = subprocess.run(
        ["samtools", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    
    if result.returncode != 0:
        self.logger.critical("PARSER: Error: samtools not found. Please install samtools and try again.")
        sys.exit(1)
    self.logger.info("PARSER: samtools was found! Proceeding with parsing")
  
  def run(self):
    """
    This function runs the parsing module for the tb_parser tool.
    """    
    self.logger.info("PARSER: Checking for samtools...")
    self.check_dependency_exists
    
    self.logger.info("PARSER: Creating initial coverage report")
    coverage = Coverage(self.logger, self.input_bam, self.output_prefix, self.coverage_regions, self.tngs, self.tngs_expert_regions)
    coverage.calculate_coverage()
    
    if self.tngs:
      self.logger.info("PARSER: Calculating the coverage for the expert rule regions")
      coverage.calculate_r_expert_rule_regions_coverage()
    
    self.logger.info("PARSER: Creating laboratorian report")
    laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix, self.tngs)
    laboratorian.create_laboratorian_report()
    
    self.logger.info("PARSER: Creating LIMS report")
    lims = LIMS(self.logger, self.input_json, self.output_prefix, self.tngs)
    lims.create_lims_report()
    
    self.logger.info("PARSER: Creating Looker report")
    looker = Looker(self.logger, self.input_json, self.output_prefix)
    looker.create_looker_report()
    
    self.logger.info("PARSER: Finalizing coverage report")
    coverage.reformat_coverage()
    
    self.logger.info("PARSER: Parsing completed")