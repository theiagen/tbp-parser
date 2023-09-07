from tbprofiler_parser.Variant import Variant

class TestVariant:
  def test_get_position_nucleotide(self):
    assert Variant.get_position(None, "c.1234A>T") == 1234
  
  def test_get_position_aminoacid(self):
    assert Variant.get_position(None, "p.Met291Ile") == 291
  