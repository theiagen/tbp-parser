import tbprofiler_parser.globals as globals

class TestGlobals:
  def test_get_position_nucleotide(self):
    assert globals.get_position("c.1234A>T") == 1234
  
  def test_get_position_aminoacid(self):
    assert globals.get_position("p.Met291Ile") == 291
  
  def test_get_position_blank(self):
    assert globals.get_position("") == 0
  
  def test_get_position_nonsensical(self):
    assert globals.get_position("asdf") == 0
