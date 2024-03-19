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
  COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, coverage_regions=COVERAGE_BED, output_prefix="test")
  COVERAGE1.calculate_coverage()
  
  INPUT_JSON = os.path.join(data_dir, "combined.json")
  
  # def test_get_lineage_and_id_bcg(self):
  #   JSON = os.path.join(self.data_dir + '/lineages', "bcg.json")
    
  #   LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
  #   lineage, ID = LOOKER1.get_lineage_and_id()
    
  #   assert (lineage, ID) == ("NA",
  #                            "M. bovis BCG")

  # def test_get_lineage_and_id_bovis(self): 
  #   JSON = os.path.join(self.data_dir + '/lineages', "bovis.json")
    
  #   LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
  #   lineage, ID = LOOKER1.get_lineage_and_id()
    
  #   assert (lineage, ID) == ("NA",
  #                            "M. bovis, not BCG")

  # def test_get_lineage_and_id_la1(self):
  #   JSON = os.path.join(self.data_dir + '/lineages', "la1.json")
    
  #   LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
  #   lineage, ID = LOOKER1.get_lineage_and_id()
    
  #   assert (lineage, ID) == ("NA",
  #                            "NA")
  
  # def test_get_lineage_and_id_lineage(self):
  #   JSON = os.path.join(self.data_dir + '/lineages', "lineage.json")
    
  #   LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
  #   lineage, ID = LOOKER1.get_lineage_and_id()
    
  #   assert (lineage, ID) == ("lineage1.1.2",
  #                            "MtBC, not M. bovis")
  
  # def test_get_lineage_and_id_nolineage(self):  
  #   JSON = os.path.join(self.data_dir + '/lineages', "nolineage.json")
    
  #   LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
  #   lineage, ID = LOOKER1.get_lineage_and_id()
    
  #   assert (lineage, ID) == ("NA", "NA")
  
  def test_create_looker_report(self):
    LOOKER1 = Looker(logger=self.LOGGER, input_json=self.INPUT_JSON, output_prefix="test")
    LOOKER1.create_looker_report()
    
    assert os.path.exists("test.looker_report.csv")
    
    os.remove("test.looker_report.csv")
    