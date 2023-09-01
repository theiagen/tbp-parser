import logging
from Coverage import Coverage
from Laboratorian import Laboratorian


class Parser:
  """
  This class runs the parsing module for the TBProfiler_Parser tool.
  """
  def __init__(self, options):
    logging.basicConfig(encoding='utf-8', level=logging.DEBUG)
    self.logger = logging.getLogger(__name__)
    self.input_json = options.input_json
    self.input_bam = options.input_bam
    self.verbose = options.verbose
    self.output_prefix = options.output_prefix
    self.min_depth = options.min_depth
    self.coverage_threshold = options.coverage_threshold
    
    if self.verbose:
      print("should be here")
      self.logger.setLevel(logging.DEBUG)
      self.logger.info("Verbose mode enabled")
    else:
      self.logger.setLevel(logging.ERROR)
  
  def run(self):
    """
    This function runs the parsing module for the TBProfiler_Parser tool.
    """    
    self.logger.info("Creating initial coverage report")
    coverage = Coverage(self.logger, self.input_bam, self.output_prefix)
    coverage.calculate_coverage(self.min_depth)
    
    self.logger.info("Creating laboratorian report")
    laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix)
    laboratorian.create_laboratorian_report(self.min_depth, self.coverage_threshold)