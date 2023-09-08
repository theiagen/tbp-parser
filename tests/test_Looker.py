import os
import logging
from tbprofiler_parser.Looker import Looker
from tbprofiler_parser.Coverage import Coverage

class TestLIMS:
  test_modules_dir = os.path.dirname(os.path.realpath(__file__))
  data_dir = os.path.join(test_modules_dir, "data")
  LOGGER = logging.getLogger(__name__)

  BAM = os.path.join(data_dir, "mtb.bam")
  COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, output_prefix="test1")
  COVERAGE1.calculate_coverage()
  
  def test_get_lineage_and_id_bcg(self):
    JSON = os.path.join(self.data_dir + '/lineages', "bcg.json")
    
    LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage, ID = LOOKER1.get_lineage_and_id()
    
    assert (lineage, ID) == ("NA",
                             "M. bovis BCG")

  def test_get_lineage_and_id_bovis(self): 
    JSON = os.path.join(self.data_dir + '/lineages', "bovis.json")
    
    LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage, ID = LOOKER1.get_lineage_and_id()
    
    assert (lineage, ID) == ("NA",
                             "M. bovis, not BCG")

  def test_get_lineage_and_id_la1(self):
    JSON = os.path.join(self.data_dir + '/lineages', "la1.json")
    
    LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage, ID = LOOKER1.get_lineage_and_id()
    
    assert (lineage, ID) == ("NA",
                             "NA")
  
  def test_get_lineage_and_id_lineage(self):
    JSON = os.path.join(self.data_dir + '/lineages', "lineage.json")
    
    LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage, ID = LOOKER1.get_lineage_and_id()
    
    assert (lineage, ID) == ("lineage1.1.2",
                             "MtBC, not M. bovis")
  
  def test_get_lienage_and_id_nolineage(self):
    JSON = os.path.join(self.data_dir + '/lineages', "nolineage.json")
    
    LOOKER1 = Looker(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage, ID = LOOKER1.get_lineage_and_id()
    
    assert (lineage, ID) == ("NA",
                             "NA")