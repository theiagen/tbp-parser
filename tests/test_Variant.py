import json
import os
import logging
from tbprofiler_parser.Variant import Variant

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, "data")

class TestVariant:
  def test_get_position_nucleotide(self):
    assert Variant.get_position(None, "c.1234A>T") == 1234
  
  def test_get_position_aminoacid(self):
    assert Variant.get_position(None, "p.Met291Ile") == 291
  
  def test_get_position_blank(self):
    assert Variant.get_position(None, "") == 0
  
  def test_get_position_nonsensical(self):
    assert Variant.get_position(None, "asdf") == 0

  def test_variant_contructor(self):
    with open(os.path.join(data_dir + "/rule1", "rrl.json"), "r") as rrl_uncertain:
      testA = json.loads(rrl_uncertain.read())
      VARIANT1 = Variant(logger=logging.getLogger(__name__), variant=testA)

    assert isinstance(VARIANT1, Variant)
    assert isinstance(VARIANT1.annotation_dictionary, dict)
    
  def test_expert_rule_rrl(self):
    with open(os.path.join(data_dir + "/rule1", "rrl.json"), "r") as rrl_uncertain:
      rrl_file = json.loads(rrl_uncertain.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in rrl_file["Rule1_Variants"]:
        rrl_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(rrl_variant.apply_expert_rules("mdl"))
        results["looker"].append(rrl_variant.apply_expert_rules("looker"))
        
    assert (results["looker"], results["mdl"]) == (["U", "U"], ["U", "S"])
  
  def test_expert_rule_rv0678(self):
    with open(os.path.join(data_dir + "/rule1", "Rv0678.json"), "r") as rv0678_uncertain:
      rv0678_file = json.loads(rv0678_uncertain.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in rv0678_file["Rule1_Variants"]:
        rv0678_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(rv0678_variant.apply_expert_rules("mdl"))
        results["looker"].append(rv0678_variant.apply_expert_rules("looker"))
      print(results)
    assert (results["looker"], results["mdl"]) == (["U", "U", "U"], ["U", "S", "U"])
      