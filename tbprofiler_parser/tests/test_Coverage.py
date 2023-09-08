import json
import os
import logging
from tbprofiler_parser.Coverage import Coverage

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, "data")

LOGGER = logging.getLogger(__name__)

BAM = os.path.join(data_dir, "mtb.bam")

COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, output_prefix="test1")

class TestCoverage:
  def test_calculate_coverage(self):
    assert isinstance(COVERAGE1, Coverage)