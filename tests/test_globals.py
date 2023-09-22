import tbp_parser.globals as globals

class TestGlobals:
  def test_get_position_nucleotide(self):
    assert globals.get_position("c.1234A>T")[0] == 1234
  
  def test_get_position_aminoacid(self):
    assert globals.get_position("p.Met291Ile")[0] == 291
  
  def test_get_position_blank(self):
    assert globals.get_position("")[0] == 0
  
  def test_get_position_nonsensical(self):
    assert globals.get_position("asdf")[0] == 0
    
  def test_get_two_positions(self):
    assert globals.get_position("c.1234_1235insA") == [1234, 1235]
