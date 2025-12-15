import os
import logging
import tbp_parser.globals
from tbp_parser.LIMS import LIMS
from tbp_parser.Coverage import Coverage
from tbp_parser.Laboratorian import Laboratorian

class TestLIMS:
    test_modules_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(test_modules_dir, "data")
    LOGGER = logging.getLogger(__name__)
    LOGGER.setLevel(logging.DEBUG)

    BAM = os.path.join(data_dir, "mtb.bam")
    COVERAGE_BED = os.path.join(data_dir, "tbdb-modified-regions-for-tests.bed")
    COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, output_prefix="test", coverage_regions=COVERAGE_BED)
    COVERAGE_DICTIONARY, AVERAGE_LOCUS_COVERAGE, LOW_DEPTH_OF_COVERAGE_LIST = COVERAGE1.get_coverage(1.0, 10)


    INPUT_JSON = os.path.join(data_dir, "combined.json")
    LABORATORIAN = Laboratorian(logger=LOGGER, input_json=INPUT_JSON, output_prefix="test")
    LABORATORIAN.create_laboratorian_report()

    def test_get_id_bcg(self):
        JSON = os.path.join(self.data_dir + '/lineages', "bcg.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id(test=True)

        assert lineage == "DNA of Mycobacterium bovis BCG detected"

    def test_get_id_bovis(self): 
        JSON = os.path.join(self.data_dir + '/lineages', "bovis.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id(test=True)

        assert lineage == "DNA of Mycobacterium bovis (not BCG) detected"

    def test_get_id_la1(self):
        JSON = os.path.join(self.data_dir + '/lineages', "la1.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id(test=True)

        assert lineage == "DNA of Mycobacterium bovis (not BCG) detected"

    def test_get_id_la1andbcg(self):
        JSON = os.path.join(self.data_dir + '/lineages', "la1andbcg.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id(test=True)

        assert lineage == "DNA of Mycobacterium bovis (not BCG) detected; DNA of Mycobacterium bovis BCG detected"

    def test_get_id_lineage(self):
        JSON = os.path.join(self.data_dir + '/lineages', "lineage.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id(test=True)

        assert lineage == "DNA of Mycobacterium tuberculosis species detected"

    def test_get_id_lineageandla1(self):
        JSON = os.path.join(self.data_dir + '/lineages', "lineageandla1.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id(test=True)

        assert lineage == "DNA of Mycobacterium bovis (not BCG) detected; DNA of Mycobacterium tuberculosis species detected"

    def test_get_id_nolineage(self):
        JSON = os.path.join(self.data_dir + '/lineages', "nolineage.json")

        LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
        lineage = LIMS1.get_id() # not setting test to check the default failing condition

        assert lineage == "DNA of Mycobacterium tuberculosis complex NOT detected"

    def test_convert_annotation(self):
        LIMS1 = LIMS(logger=self.LOGGER, input_json=self.INPUT_JSON, output_prefix="test")

        message1 = LIMS1.convert_annotation("R", "rifampicin")
        message2 = LIMS1.convert_annotation("S", "isoniazid")
        message3 = LIMS1.convert_annotation("R-Interim", "ethambutol")
        message4 = LIMS1.convert_annotation("Insufficient Coverage", "kanamycin")

        assert (message1, message2, message3, message4) == ("Mutation(s) associated with resistance to rifampicin detected", 
                                                  "No mutations associated with resistance to isoniazid detected", 
                                                  "The detected mutation(s) have uncertain significance. Resistance to ethambutol cannot be ruled out",
                                                  "Pending Retest")

    def test_create_lims_report(self):   
        LIMS1 = LIMS(logger=self.LOGGER, input_json=self.INPUT_JSON, output_prefix="test")
        LIMS1.create_lims_report()

        assert os.path.exists("test.lims_report.csv")

        os.remove("test.lims_report.csv")