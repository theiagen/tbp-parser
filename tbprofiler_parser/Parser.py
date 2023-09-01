import logging
import os
import sys
import time
from tbprofiler_parser.Coverage import Coverage
from tbprofiler_parser.Laboratorian import Laboratorian


class Parser:
  """
  This class runs the parsing module for the TBProfiler_Parser tool.
  """
  def __init__(self, options):
    self.logger = logging.getLogger(__name__)
    self.logger.info("Initializing Parser")
    self.input_json = options.input_json
    self.input_bam = options.input_bam
    self.verbose = options.verbose
    self.output_prefix = options.output_prefix
    self.min_depth = options.min_depth
    self.coverage_threshold = options.coverage_threshold
    
    if self.verbose:
      self.logger.info("Verbose mode enabled")
      self.logger.setLevel(logging.DEBUG)
    else:
      self.logger.setLevel(logging.ERROR)
  
  def run(self):
    """
    This function runs the parsing module for the TBProfiler_Parser tool.
    """
    self.logger.info("Creating initial coverage report")
    coverage = Coverage(self.logger, self.input_bam, self.output_prefix)
    coverage.calculate_coverage()
    
    self.logger.info("Creating laboratorian report")
    laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix)
    