import os
import logging
from tbp_parser.Looker import Looker
from tbp_parser.Coverage import Coverage

class TestLIMS:
    test_modules_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(test_modules_dir, "data")
    LOGGER = logging.getLogger(__name__)

    BAM = os.path.join(data_dir, "mtb.bam")
    COVERAGE_BED = os.path.join(data_dir, "tbdb-modified-regions-for-tests.bed")
    
    COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, OUTPUT_PREFIX="test", coverage_regions=COVERAGE_BED)
    COVERAGE_DICTIONARY, AVERAGE_LOCUS_COVERAGE, LOW_DEPTH_OF_COVERAGE_LIST = COVERAGE1.get_coverage(1.0, 10)


    def test_create_looker_report(self):
        LOOKER1 = Looker(logger=self.LOGGER, OUTPUT_PREFIX="test")
        LOOKER1.create_looker_report()

        assert os.path.exists("test.looker_report.csv")

        os.remove("test.looker_report.csv")
    