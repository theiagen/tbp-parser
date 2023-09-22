from tbp_parser.Row import Row
import logging

class TestRow:  
  def test_rank_annotation(self):
    r = Row(logging.getLogger(__name__), None, "Assoc w R", None, "test").rank_annotation()
    r_interim = Row(logging.getLogger(__name__), None, "Assoc w R - interim", None, "test").rank_annotation()
    u = Row(logging.getLogger(__name__), None, "Uncertain significance", None, "test").rank_annotation()
    s_interim = Row(logging.getLogger(__name__), None, "Not assoc w R - Interim", None, "test").rank_annotation()
    
    assert (r, r_interim, u, s_interim) == (4, 3, 2, 1)
  
  def test_annotation_to_LIMS(self):
    r = Row(logging.getLogger(__name__), None, "Assoc w R", "test", "test").annotation_to_LIMS()
    r_interim = Row(logging.getLogger(__name__), None, "Assoc w R - interim", "test", "test").annotation_to_LIMS()
    u = Row(logging.getLogger(__name__), None, "Uncertain significance", "test", "test").annotation_to_LIMS()
    s_interim = Row(logging.getLogger(__name__), None, "Not assoc w R - Interim", "test", "test").annotation_to_LIMS()
    
    assert(r, r_interim, u, s_interim) == ("Mutation(s) associated with resistance to test detected", 
                                          "The detected mutation(s) have uncertain significance. Resistance to test cannot be ruled out", 
                                          "The detected mutation(s) have uncertain significance. Resistance to test cannot be ruled out", 
                                          "No mutations associated with resistance to test detected")
    
  def test_remove_no_expert(self):
    a = Row(logging.getLogger(__name__), None, "Assoc w R", "test", "test")
    a.looker_interpretation = "Unoexpert"
    a.mdl_interpretation = "Snoexpert"
    
    a.remove_no_expert()
    
    assert (a.looker_interpretation, a.mdl_interpretation, a.rationale, a.confidence) == ("U", "S", "No WHO annotation or expert rule", "No WHO annotation")
    
  def test_complete_row(self):
    a = Row(logging.getLogger(__name__), None, "Assoc w R", "test", "katG")
    a.complete_row()
    assert a.rationale == "WHO classification"