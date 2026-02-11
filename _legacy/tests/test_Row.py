from tbp_parser.Row import Row
import logging

class TestRow:  
    def test_rank_annotation(self):
        r = Row(logging.getLogger(__name__), None, "Assoc w R", None, "test").rank_annotation()
        r_interim = Row(logging.getLogger(__name__), None, "Assoc w R - interim", None, "test").rank_annotation()
        u = Row(logging.getLogger(__name__), None, "Uncertain significance", None, "test").rank_annotation()
        s_interim = Row(logging.getLogger(__name__), None, "Not assoc w R - Interim", None, "test").rank_annotation()

        assert (r, r_interim, u, s_interim) == (4, 3, 2, 1)

    def test_describe_rationale_rule324(self):
        a = Row(logging.getLogger(__name__), None, "Assoc w R", "test", "test")
        a.looker_interpretation = "Urule3.2.4"
        a.mdl_interpretation = "Srule3.2.4"
        a.rationale = ""
        a.confidence = ""

        a.describe_rationale()

        assert (a.looker_interpretation, a.mdl_interpretation, a.rationale, a.confidence) == ("U", "S", "No WHO annotation or expert rule", "No WHO annotation")

    def test_describe_rationale_whov2(self):
        a = Row(logging.getLogger(__name__), None, "", "test", "test")
        a.looker_interpretation = "R", "whov2"
        a.mdl_interpretation = "Uwhov2"
        a.rationale = ""
        a.confidence = ""

        a.describe_rationale()

        assert (a.looker_interpretation, a.mdl_interpretation, a.rationale, a.confidence) == ("R", "U", "Mutation in proximal promoter region", "No WHO annotation")

    def test_complete_row(self):
        a = Row(logging.getLogger(__name__), None, "Assoc w R", "test", "katG")
        a.complete_row()
        assert a.rationale == "WHO classification"