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

  def test_is_within_range(self):
    assert globals.is_within_range([1234], [1233, 1235]) == True
    
  def test_is_not_within_range(self):
    assert globals.is_within_range([1236], [1234, 1235]) == False
    
  def test_is_within_range_deletion(self):
    assert globals.is_within_range([1234, 1236], [1233, 1235]) == True
    
  def test_is_not_within_range_deletion(self):
    assert globals.is_within_range([1234, 1236], [1237, 1238]) == False
    
  def test_is_within_two_ranges(self):
    assert globals.is_within_range([1234, 1236], [[1233, 1235], [1237, 1238]]) == True
  
  def test_is_not_within_two_ranges(self):
    assert globals.is_within_range([1234, 1236], [[1232, 1233], [1237, 1238]]) == False