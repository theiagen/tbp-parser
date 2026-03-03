import pytest
from variant import VariantInterpreter, InterpretationResult


@pytest.fixture
def interpreter() -> VariantInterpreter:
    """Create a VariantInterpreter"""
    return VariantInterpreter()

class TestDetermineInterpretation:
    """Tests for the new determine_interpretation method."""

    def test_who_annotated_variant(self, interpreter, make_variant):
        v = make_variant(confidence="Assoc w R")
        interpreter.determine_interpretation([v])

        assert v.mdl_interpretation == "R"
        assert v.looker_interpretation == "R"
        assert "WHO classification" in v.rationale

    def test_interim_who_annotation(self, interpreter, make_variant):
        v = make_variant(confidence="Assoc w R - Interim")
        interpreter.determine_interpretation([v])

        assert v.mdl_interpretation == "R"
        assert v.looker_interpretation == "R-Interim"

    def test_no_who_annotation_applies_expert_rules(self, interpreter, make_variant):
        v = make_variant(
            gene_name="mmpR5", gene_id="Rv0678",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Ala",
        )
        interpreter.determine_interpretation([v])

        assert v.mdl_interpretation == "U"
        assert v.rationale is not None
        assert "1.2" in v.rationale  # Should use rule 1.2 for CDC genes


class TestInterpret:
    """Tests for the new interpret() method which returns InterpretationResult."""

    # =========================================================================
    # Section 1: CDC Expert Rule Genes (mmpR5, atpE, pepQ, mmpL5, mmpS5, rrl, rplC)
    # =========================================================================

    def test_rule1_1_cdc_gene_with_who_confidence(self, interpreter, make_variant):
        """CDC gene (mmpR5) with WHO confidence → rule 1.1"""
        v = make_variant(
            gene_name="mmpR5", gene_id="Rv0678",
            confidence="Assoc w R",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Ala",
        )
        result = interpreter.interpret(v)

        assert isinstance(result, InterpretationResult)
        assert result.rule_id == "1.1"
        assert result.looker_interpretation == "R"
        assert result.mdl_interpretation == "R"

    def test_rule1_2_novel_drug_target_nonsynonymous(self, interpreter, make_variant):
        """mmpR5 nonsynonymous mutation without WHO → U, rule 1.2"""
        v = make_variant(
            gene_name="mmpR5", gene_id="Rv0678",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Ala",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "U"
        assert result.mdl_interpretation == "U"
        assert result.rule_id == "1.2"

    def test_rule1_2_novel_drug_target_synonymous(self, interpreter, make_variant):
        """mmpR5 synonymous mutation without WHO → S, rule 1.2"""
        v = make_variant(
            gene_name="mmpR5", gene_id="Rv0678",
            confidence="No WHO annotation",
            type="synonymous_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Thr",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "1.2"

    def test_rule1_2_upstream_gene_variant(self, interpreter, make_variant):
        """mmpR5 upstream gene variant without WHO → S, rule 1.2"""
        v = make_variant(
            gene_name="mmpR5", gene_id="Rv0678",
            confidence="No WHO annotation",
            type="upstream_gene_variant",
            nucleotide_change="c.-100A>G",
            protein_change="",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "1.2"

    def test_rule1_2_rv0678_also_recognized(self, interpreter, make_variant):
        """Rv0678 (mmpR5 locus tag) is also recognized as CDC gene"""
        v = make_variant(
            gene_name="Rv0678", gene_id="Rv0678",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Ala",
        )
        result = interpreter.interpret(v)

        assert result.rule_id == "1.2"  # Should use CDC rules, not section 3

    # =========================================================================
    # Section 2: WHO Expert Rule Genes (katG, pncA, ethA, gid, rpoB)
    # =========================================================================

    def test_rule2_1_who_gene_with_who_confidence(self, interpreter, make_variant):
        """WHO gene (katG) with WHO confidence → rule 2.1"""
        v = make_variant(
            gene_name="katG", gene_id="Rv1908c",
            confidence="Assoc w R",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Ala",
        )
        result = interpreter.interpret(v)

        assert result.rule_id == "2.1"
        assert result.looker_interpretation == "R"

    def test_rule2_2_1_1_katg_lof(self, interpreter, make_variant):
        """katG loss-of-function (deletion) within 30 nt → R, rule 2.2.1.1"""
        v = make_variant(
            gene_name="katG", gene_id="Rv1908c",
            confidence="No WHO annotation",
            type="frameshift_variant",
            nucleotide_change="c.1_10del",
            protein_change="p.Met1fs",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "R"
        assert result.mdl_interpretation == "R"
        assert result.rule_id == "2.2.1.1"

    def test_rule2_2_1_1_katg_synonymous(self, interpreter, make_variant):
        """katG synonymous (not LOF) → S, rule 2.2.1.1"""
        v = make_variant(
            gene_name="katG", gene_id="Rv1908c",
            confidence="No WHO annotation",
            type="synonymous_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Thr",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "2.2.1.1"

    def test_rule2_2_2_1_rpob_rrdr_nonsynonymous(self, interpreter, make_variant):
        """rpoB RRDR nonsynonymous → R, rule 2.2.2.1"""
        v = make_variant(
            gene_name="rpoB", gene_id="Rv0667",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.1349C>T",
            protein_change="p.Ser450Leu",  # codon 450 is in RRDR [426,452]
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "R"
        assert result.mdl_interpretation == "R"
        assert result.rule_id == "2.2.2.1"

    def test_rule2_2_2_1_rpob_rrdr_synonymous(self, interpreter, make_variant):
        """rpoB RRDR synonymous → S, rule 2.2.2.1"""
        v = make_variant(
            gene_name="rpoB", gene_id="Rv0667",
            confidence="No WHO annotation",
            type="synonymous_variant",
            nucleotide_change="c.1350C>T",
            protein_change="p.Ser450Ser",  # codon 450 is in RRDR
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "2.2.2.1"

    def test_rule2_2_2_2_rpob_non_rrdr_nonsynonymous(self, interpreter, make_variant):
        """rpoB outside RRDR, nonsynonymous → U, rule 2.2.2.2"""
        v = make_variant(
            gene_name="rpoB", gene_id="Rv0667",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Ala",  # codon 34 is outside RRDR
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "U"
        assert result.mdl_interpretation == "U"
        assert result.rule_id == "2.2.2.2"

    def test_rule2_2_2_2_rpob_non_rrdr_synonymous(self, interpreter, make_variant):
        """rpoB outside RRDR, synonymous → S, rule 2.2.2.2"""
        v = make_variant(
            gene_name="rpoB", gene_id="Rv0667",
            confidence="No WHO annotation",
            type="synonymous_variant",
            nucleotide_change="c.100A>G",
            protein_change="p.Thr34Thr",  # codon 34 is outside RRDR
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "2.2.2.2"

    # =========================================================================
    # Section 3: Other Genes (rrs, gyrA, gyrB, and default)
    # =========================================================================

    def test_rule3_1_other_gene_with_who_confidence(self, interpreter, make_variant):
        """Other gene (embB) with WHO confidence → rule 3.1"""
        v = make_variant(
            gene_name="embB", gene_id="Rv3795",
            confidence="Assoc w R",
            type="missense_variant",
            nucleotide_change="c.916A>G",
            protein_change="p.Met306Val",
        )
        result = interpreter.interpret(v)

        assert result.rule_id == "3.1"
        assert result.looker_interpretation == "R"

    def test_rule3_2_1_rrs_special_position(self, interpreter, make_variant):
        """rrs at special position (1401) → U, rule 3.2.1"""
        v = make_variant(
            gene_name="rrs", gene_id="EBG00000313325",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.1401A>G",
            protein_change="",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "U"
        assert result.mdl_interpretation == "U"
        assert result.rule_id == "3.2.1"

    def test_rule3_2_1_rrs_non_special_position(self, interpreter, make_variant):
        """rrs at non-special position → S, rule 3.2.1"""
        v = make_variant(
            gene_name="rrs", gene_id="EBG00000313325",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.100A>G",
            protein_change="",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "3.2.1"

    def test_rule3_2_2_gyra_qrdr_nonsynonymous(self, interpreter, make_variant):
        """gyrA QRDR nonsynonymous → U, rule 3.2.2"""
        v = make_variant(
            gene_name="gyrA", gene_id="Rv0006",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.271G>C",
            protein_change="p.Asp91His",  # codon 91 is in QRDR [88,94]
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "U"
        assert result.mdl_interpretation == "U"
        assert result.rule_id == "3.2.2"

    def test_rule3_2_3_gyrb_qrdr_nonsynonymous(self, interpreter, make_variant):
        """gyrB QRDR nonsynonymous → U, rule 3.2.3"""
        v = make_variant(
            gene_name="gyrB", gene_id="Rv0005",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.1400A>G",
            protein_change="p.Asp467Gly",  # codon 467 is in QRDR [446,507]
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "U"
        assert result.mdl_interpretation == "U"
        assert result.rule_id == "3.2.3"

    def test_rule3_2_4_other_gene_nonsynonymous(self, interpreter, make_variant):
        """Gene not in special lists, nonsynonymous → U, rule 3.2.4"""
        v = make_variant(
            gene_name="embB", gene_id="Rv3795",
            confidence="No WHO annotation",
            type="missense_variant",
            nucleotide_change="c.916A>G",
            protein_change="p.Met306Val",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "U"
        assert result.mdl_interpretation == "U"
        assert result.rule_id == "3.2.4"

    def test_rule3_2_4_other_gene_synonymous(self, interpreter, make_variant):
        """Gene not in special lists, synonymous → S, rule 3.2.4"""
        v = make_variant(
            gene_name="embB", gene_id="Rv3795",
            confidence="No WHO annotation",
            type="synonymous_variant",
            nucleotide_change="c.900C>T",
            protein_change="p.Ala300Ala",
        )
        result = interpreter.interpret(v)

        assert result.looker_interpretation == "S"
        assert result.mdl_interpretation == "S"
        assert result.rule_id == "3.2.4"


class TestInterpretationResult:
    """Tests for the InterpretationResult dataclass."""

    def test_interpretation_result_fields(self):
        """Verify InterpretationResult has all expected fields."""
        result = InterpretationResult(
            confidence="Assoc w R",
            looker_interpretation="R",
            mdl_interpretation="R",
            rule_id="1.1",
            rationale="WHO confidence available"
        )

        assert result.confidence == "Assoc w R"
        assert result.looker_interpretation == "R"
        assert result.mdl_interpretation == "R"
        assert result.rule_id == "1.1"
        assert result.rationale == "WHO confidence available"
