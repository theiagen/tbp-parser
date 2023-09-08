import json
import os
import logging
from tbprofiler_parser.Variant import Variant

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, "data")

LOGGER = logging.getLogger(__name__)

with open(os.path.join(data_dir, "TestA.json"), "r") as f:
  testA = json.loads(f.read())
  VARIANT1 = Variant(logger=LOGGER, variant=testA)


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
    assert isinstance(VARIANT1, Variant)
    assert isinstance(VARIANT1.annotation_dictionary, dict)