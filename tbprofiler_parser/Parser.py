import logging
import globals
import sys
from Coverage import Coverage
from Laboratorian import Laboratorian
from Looker import Looker
from LIMS import LIMS


class Parser:
  """
  This class runs the parsing module for the TBProfiler_Parser tool.
  """
  def __init__(self, options):
    logging.basicConfig(encoding='utf-8', level=logging.ERROR, stream=sys.stderr)
    self.logger = logging.getLogger(__name__)
    self.input_json = options.input_json
    self.input_bam = options.input_bam
    self.verbose = options.verbose
    self.debug = options.debug
    self.output_prefix = options.output_prefix
    globals.MIN_DEPTH = options.min_depth
    globals.COVERAGE_THRESHOLD = options.coverage_threshold
    globals.SEQUENCING_METHOD = options.sequencing_method
    globals.OPERATOR = options.operator
    
    if self.verbose:
      self.logger.setLevel(logging.INFO)
      self.logger.info("Verbose mode enabled")
    else:
      self.logger.setLevel(logging.ERROR)
      
    if self.debug:
      self.logger.setLevel(logging.DEBUG)
      self.logger.debug("Debug mode enabled")
  
  def run(self):
    """
    This function runs the parsing module for the TBProfiler_Parser tool.
    """    
    self.logger.info("Creating initial coverage report")
    coverage = Coverage(self.logger, self.input_bam, self.output_prefix)
    coverage.calculate_coverage()
    
    self.logger.info("Creating laboratorian report")
    laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix)
    laboratorian.create_laboratorian_report()
    
    self.logger.info("Creating Looker report")
    looker = Looker(self.logger, self.input_json, self.output_prefix)
    looker.create_looker_report()
    
    self.logger.info("Creating LIMS report")
    lims = LIMS(self.logger, self.input_json, self.output_prefix)
    lims.create_lims_report()
    
    self.logger.info("Finalizing coverage report")
    coverage.reformat_coverage()
    
    self.logger.info("Parsing completed")