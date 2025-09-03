from tbp_parser.Laboratorian import Laboratorian
from tbp_parser.Coverage import Coverage
import logging
import os
import json

class TestLaboratorian:
  test_modules_dir = os.path.dirname(os.path.realpath(__file__))
  data_dir = os.path.join(test_modules_dir, "data")
  LOGGER = logging.getLogger(__name__)
  JSON = os.path.join(data_dir, "combined.json")

  BAM = os.path.join(data_dir, "mtb.bam")
  COVERAGE_BED = os.path.join(data_dir, "tbdb-modified-regions-for-tests.bed")
  COVERAGE1 = Coverage(logger=LOGGER, input_bam=BAM, coverage_regions=COVERAGE_BED, output_prefix="test", tngs_expert_regions=None)
  COVERAGE1.calculate_coverage()
  
  def test_iterate_section(self):
    with open(self.JSON, "r") as input:
      input_section = json.loads(input.read())
      
      LAB1 = Laboratorian(logger=self.LOGGER, input_json=self.JSON, output_prefix="test")
      
      row_list = LAB1.iterate_section(input_section["Rule1_Variants"], [])
      
    assert len(row_list) == 8
      
  def test_create_lab_report(self):
    LAB1 = Laboratorian(logger=self.LOGGER, input_json=self.JSON, output_prefix="test")
    LAB1.create_laboratorian_report()    
    
    assert os.path.exists("test.laboratorian_report.csv")
    
    os.remove("test.laboratorian_report.csv")