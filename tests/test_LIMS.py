import os
import logging
from tbp_parser.LIMS import LIMS
from tbp_parser.Coverage import Coverage
from tbp_parser.Laboratorian import Laboratorian

class TestLIMS:
  test_modules_dir = os.path.dirname(os.path.realpath(__file__))
  data_dir = os.path.join(test_modules_dir, "data")
  LOGGER = logging.getLogger(__name__)

  BAM = os.path.join(data_dir, "mtb.bam")
  COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, output_prefix="test")
  COVERAGE1.calculate_coverage()
  
  INPUT_JSON = os.path.join(data_dir, "combined.json")
  LABORATORIAN = Laboratorian(logger=LOGGER, input_json=INPUT_JSON, output_prefix="test")
  LABORATORIAN.create_laboratorian_report()
  
  def test_get_lineage_bcg(self):
    JSON = os.path.join(self.data_dir + '/lineages', "bcg.json")
    
    LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage = LIMS1.get_lineage()
    
    assert lineage == "DNA of Mycobacterium bovis BCG detected"

  def test_get_lineage_bovis(self): 
    JSON = os.path.join(self.data_dir + '/lineages', "bovis.json")
    
    LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage = LIMS1.get_lineage()
    
    assert lineage == "DNA of non-BCG Mycobacterium bovis detected"

  def test_get_lineage_la1(self):
    JSON = os.path.join(self.data_dir + '/lineages', "la1.json")
    
    LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage = LIMS1.get_lineage()
    
    assert lineage == "DNA of M. tuberculosis complex detected (M. bovis)"
  
  def test_get_lineage_lineage(self):
    JSON = os.path.join(self.data_dir + '/lineages', "lineage.json")
    
    LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage = LIMS1.get_lineage()
    
    assert lineage == "DNA of Mycobacterium tuberculosis species detected"
  
  def test_get_lineage_nolineage(self):
    JSON = os.path.join(self.data_dir + '/lineages', "nolineage.json")
    
    LIMS1 = LIMS(logger=self.LOGGER, input_json=JSON, output_prefix="test")
    lineage = LIMS1.get_lineage()
    
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