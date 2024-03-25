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
    globals.COVERAGE_THRESHOLD = options.coverage_threshold
    globals.SEQUENCING_METHOD = options.sequencing_method
    globals.RRS_FREQUENCY = options.rrs_frequency
    globals.RRL_FREQUENCY = options.rrl_frequency
    globals.OPERATOR = options.operator
    
    if self.tngs:
      self.logger.info("Deeplex + CDPH modified protocol flag detected; adjusting outputs to reflect this")
      if (self.coverage_regions == "../data/tbdb-modified-regions.bed"):
        self.logger.debug("Changing default coverage regions to ../data/tngs-reportable-regions.bed")
        self.coverage_regions = "../data/tngs-reportable-regions.bed"
      
      self.logger.debug("Altering the ANTIMICROBIAL_CODE_TO_GENES dictionary to include only tNGS entries")
      globals.ANTIMICROBIAL_CODE_TO_GENES = globals.ANTIMICROBIAL_CODE_TO_GENES_tNGS
      
      self.logger.debug("Altering the GENES_FOR_LIMS list to include only tNGS genes")
      globals.GENES_FOR_LIMS = globals.GENES_FOR_LIMS_tNGS
      
    else:
      self.logger.debug("Setting the ANTIMICROBIAL_CODE_TO_GENES dictionary to include all WGS entries")
      globals.ANTIMICROBIAL_CODE_TO_GENES = globals.ANTIMICROBIAL_CODE_TO_GENES_WGS
    
      self.logger.debug("Setting the GENES_FOR_LIMS list to include all WGS genes")
      globals.GENES_FOR_LIMS = globals.GENES_FOR_LIMS_WGS
    
    if self.verbose:
      self.logger.setLevel(logging.INFO)
      self.logger.info("Verbose mode enabled")
    else:
      self.logger.setLevel(logging.ERROR)
      
    if self.debug:
      self.logger.setLevel(logging.DEBUG)
      self.logger.debug("Debug mode enabled")
  
  def check_dependency_exists(self):
    result = subprocess.run(
        ["samtools", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    
    if result.returncode != 0:
        self.logger.critical("Error: samtools not found. Please install samtools and try again.")
        sys.exit(1)
    self.logger.info("samtools was found! Proceeding with parsing")
  
  def run(self):
    """
    This function runs the parsing module for the tb_parser tool.
    """    
    self.logger.info("Checking for samtools...")
    self.check_dependency_exists
    
    self.logger.info("Creating initial coverage report")
    coverage = Coverage(self.logger, self.input_bam, self.output_prefix, self.coverage_regions, self.tngs, self.tngs_expert_regions)
    coverage.calculate_coverage()
    
    if self.tngs:
      self.logger.info("Calculating the coverage for the expert rule regions")
      coverage.calculate_r_expert_rule_regions_coverage()
    
    self.logger.info("Creating laboratorian report")
    laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix, self.tngs)
    laboratorian.create_laboratorian_report()
    
    self.logger.info("Creating LIMS report")
    lims = LIMS(self.logger, self.input_json, self.output_prefix, self.tngs)
    lims.create_lims_report()
    
    self.logger.info("Creating Looker report")
    looker = Looker(self.logger, self.input_json, self.output_prefix)
    looker.create_looker_report()
    
    self.logger.info("Finalizing coverage report")
    coverage.reformat_coverage()
    
    self.logger.info("Parsing completed")