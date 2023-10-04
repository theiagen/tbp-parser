import json
import os
import logging
from tbp_parser.Variant import Variant

class TestVariant:
  test_modules_dir = os.path.dirname(os.path.realpath(__file__))
  data_dir = os.path.join(test_modules_dir, "data")

  def test_variant_contructor(self):
    with open(os.path.join(self.data_dir + "/rule1", "rrl.json"), "r") as rrl_uncertain:
      testA = json.loads(rrl_uncertain.read())
      VARIANT1 = Variant(logger=logging.getLogger(__name__), variant=testA)

    assert isinstance(VARIANT1, Variant)
    assert isinstance(VARIANT1.annotation_dictionary, dict)
    
  def test_expert_rule_rrl(self):
    with open(os.path.join(self.data_dir + "/rule1", "rrl.json"), "r") as rrl:
      rrl_file = json.loads(rrl.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in rrl_file["Rule1_Variants"]:
        rrl_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(rrl_variant.apply_expert_rules("mdl"))
        results["looker"].append(rrl_variant.apply_expert_rules("looker"))
        
    assert (results["looker"], results["mdl"]) == (["Urule1.2", "Urule1.2"], ["Urule1.2", "Srule1.2"])
  
  def test_expert_rule_rv0678(self):
    with open(os.path.join(self.data_dir + "/rule1", "Rv0678.json"), "r") as rv0678:
      rv0678_file = json.loads(rv0678.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in rv0678_file["Rule1_Variants"]:
        rv0678_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(rv0678_variant.apply_expert_rules("mdl"))
        results["looker"].append(rv0678_variant.apply_expert_rules("looker"))
    
    assert (results["looker"], results["mdl"]) == (["Urule1.2", "Urule1.2", "Urule1.2"], ["Urule1.2", "Srule1.2", "Urule1.2"])
  
  def test_expert_rule_katg(self):
    with open(os.path.join(self.data_dir + "/rule2", "katg.json"), "r") as katg:
      katg_file = json.loads(katg.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in katg_file["Rule2_Variants"]:
        katg_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(katg_variant.apply_expert_rules("mdl"))
        results["looker"].append(katg_variant.apply_expert_rules("looker"))

    assert (results["looker"], results["mdl"]) == (["Urule2.2.1", "Urule2.2.1"], ["Urule2.2.1", "Srule2.2.1"])
  
  def test_expert_rule_rpob(self):
    with open(os.path.join(self.data_dir + "/rule2", "rpob.json"), "r") as rpob:
      rpob_file = json.loads(rpob.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in rpob_file["Rule2_Variants"]:
        rpob_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(rpob_variant.apply_expert_rules("mdl"))
        results["looker"].append(rpob_variant.apply_expert_rules("looker"))
    
    assert (results["looker"], results["mdl"]) == (["Urule2.2.2.2", "Urule2.2.2.2", "Rrule2.2.2.1"], ["Srule2.2.2.2", "Srule2.2.2.2", "Rrule2.2.2.1"])
  
  def test_expert_rule_fabg1(self):
    with open(os.path.join(self.data_dir + "/rule3", "fabg1.json"), "r") as fabg1:
      fabg1_file = json.loads(fabg1.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in fabg1_file["Rule3_Variants"]:
        fabg1_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(fabg1_variant.apply_expert_rules("mdl"))
        results["looker"].append(fabg1_variant.apply_expert_rules("looker"))
    
    assert (results["looker"], results["mdl"]) == (["Urule3.2.4"], ["Srule3.2.4"])
  
  def test_expert_rule_rrs(self):
    with open(os.path.join(self.data_dir + "/rule3", "rrs.json"), "r") as rrs:
      rrs_file = json.loads(rrs.read())
      results= {}
      results["mdl"] = []
      results["looker"] = []
      for variant in rrs_file["Rule3_Variants"]:
        rrs_variant = Variant(logger=logging.getLogger(__name__), variant=variant)
        results["mdl"].append(rrs_variant.apply_expert_rules("mdl"))
        results["looker"].append(rrs_variant.apply_expert_rules("looker"))
    
    assert (results["looker"], results["mdl"]) == (["Urule3.2.1", "Urule3.2.1"], ["Urule3.2.1", "Srule3.2.1"])