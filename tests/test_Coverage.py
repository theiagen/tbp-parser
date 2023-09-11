import os
import logging
from tbprofiler_parser.Coverage import Coverage

class TestCoverage:
  test_modules_dir = os.path.dirname(os.path.realpath(__file__))
  data_dir = os.path.join(test_modules_dir, "data")

  def test_calculate_coverage_constructor(self):
    
    LOGGER = logging.getLogger(__name__)
    BAM = os.path.join(self.data_dir, "mtb.bam")
    
    COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, output_prefix="test")

    assert isinstance(COVERAGE1, Coverage)