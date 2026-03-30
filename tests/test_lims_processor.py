import pytest
from tbp_parser.LIMS import LIMSProcessor, LIMSRecord, LIMSGeneCode, parse_lims_yml_file
from tbp_parser.Coverage.coverage_data import ERRCoverage

@pytest.fixture
def processor():
    return LIMSProcessor()

@pytest.fixture
def make_lims_gene_code():
    def _make(gene_code="M_DST_D02_rpoB"):
        return LIMSGeneCode(gene_code=gene_code)
    return _make

@pytest.fixture
def make_lims_record():
    def _make(drug="rifampicin", drug_code="M_DST_D02", gene_codes=None):
        if gene_codes is None:
            gene_codes = {"rpoB": LIMSGeneCode(gene_code="M_DST_D02_rpoB")}
        return LIMSRecord(drug=drug, drug_code=drug_code, gene_codes=gene_codes)
    return _make

class TestParseLimsYmlFile:
    def test_parse_lims_yml_num_records(self, mock_config):
        """Test that the default LIMS YAML file is parsed correctly"""
        result = parse_lims_yml_file(mock_config.lims_report_format_yml)
        assert len(result) == 18
        assert all(isinstance(r, LIMSRecord) for r in result)

class TestResistanceRanking:
    @pytest.mark.parametrize("higher,lower", [
        ("R", "Insufficient Coverage"),
        ("Insufficient Coverage", "U"),
        ("U", "S"),
        ("S", "WT"),
        ("WT", "NA"),
    ])
    def test_ranking_order(self, higher, lower):
        ranking = LIMSProcessor.RESISTANCE_RANKING
        assert ranking[higher] > ranking[lower]

class TestSetGeneCodeMaxMdl:
    def test_r_wins_over_s(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v1 = make_variant(mdl_interpretation="R")
        v2 = make_variant(mdl_interpretation="S", nucleotide_change="c.100A>G", protein_change="p.Thr34Ala")
        processor.set_gene_code_max_mdl(gc, [v1, v2])
        assert gc.max_mdl_interpretation == "R"
        assert v1 in gc.max_mdl_variants

    def test_u_wins_over_s(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v1 = make_variant(mdl_interpretation="U")
        v2 = make_variant(mdl_interpretation="S", nucleotide_change="c.100A>G", protein_change="p.Thr34Ala")
        processor.set_gene_code_max_mdl(gc, [v1, v2])
        assert gc.max_mdl_interpretation == "U"

    def test_all_u_returns_u(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v1 = make_variant(mdl_interpretation="U")
        v2 = make_variant(mdl_interpretation="U", nucleotide_change="c.100A>G", protein_change="p.Thr34Ala")
        processor.set_gene_code_max_mdl(gc, [v1, v2])
        assert gc.max_mdl_interpretation == "U"

    def test_empty_variants_returns_na(self, processor, make_lims_gene_code):
        gc = make_lims_gene_code()
        processor.set_gene_code_max_mdl(gc, [])
        assert gc.max_mdl_interpretation == "NA"

    def test_no_interpretations_returns_na(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant()
        v.mdl_interpretation = None
        processor.set_gene_code_max_mdl(gc, [v])
        assert gc.max_mdl_interpretation == "NA"

    def test_wt_only_returns_wt(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="WT")
        processor.set_gene_code_max_mdl(gc, [v])
        assert gc.max_mdl_interpretation == "WT"

    def test_r_includes_u_variants_in_reportable(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v_r = make_variant(mdl_interpretation="R")
        v_u = make_variant(mdl_interpretation="U", nucleotide_change="c.100A>G", protein_change="p.Thr34Ala")
        v_s = make_variant(mdl_interpretation="S", nucleotide_change="c.200C>T", protein_change="p.Ala67Val")
        processor.set_gene_code_max_mdl(gc, [v_r, v_u, v_s])
        assert gc.max_mdl_interpretation == "R"
        assert gc.max_mdl_variants == [v_r]
        assert v_r in gc.max_mdl_reportable_variants
        assert v_u in gc.max_mdl_reportable_variants
        assert v_s not in gc.max_mdl_reportable_variants

    def test_r_without_u_only_includes_r(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v_r = make_variant(mdl_interpretation="R")
        v_s = make_variant(mdl_interpretation="S", nucleotide_change="c.100A>G", protein_change="p.Thr34Ala")
        processor.set_gene_code_max_mdl(gc, [v_r, v_s])
        assert gc.max_mdl_interpretation == "R"
        assert gc.max_mdl_variants == [v_r]
        assert v_r in gc.max_mdl_reportable_variants
        assert v_s not in gc.max_mdl_reportable_variants

    def test_synonymous_rpob_rrdr_in_reportable_not_max(self, processor, make_lims_gene_code, make_variant):
        """Synonymous rpoB RRDR variants appear in max_mdl_reportable_variants but not max_mdl_variants when their interpretation differs from max."""
        gc = make_lims_gene_code()
        v_r = make_variant(mdl_interpretation="R", protein_change="p.Ser450Leu")
        # synonymous rpoB outside RRDR region (pos 400, outside 426-452)
        v_syn_outside = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Ala400Ala", nucleotide_change="c.1200A>G",
            type="synonymous_variant",
        )
        # synonymous rpoB inside RRDR region (pos 433, within 426-452)
        v_syn_rrdr = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Phe433Phe", nucleotide_change="c.1299C>T",
            type="synonymous_variant",
        )
        processor.set_gene_code_max_mdl(gc, [v_r, v_syn_outside, v_syn_rrdr])
        assert gc.max_mdl_interpretation == "R"
        assert gc.max_mdl_variants == [v_r]
        assert v_r in gc.max_mdl_reportable_variants
        assert v_syn_rrdr in gc.max_mdl_reportable_variants
        assert v_syn_outside not in gc.max_mdl_reportable_variants


class TestResolveGeneTarget:
    def test_r_shows_protein_change(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="R", protein_change="p.Ser450Leu")
        gc.max_mdl_interpretation = "R"
        gc.max_mdl_variants = [v]
        gc.max_mdl_reportable_variants = [v]
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "p.Ser450Leu"

    def test_u_shows_protein_change(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="U", protein_change="p.Thr34Ala")
        gc.max_mdl_interpretation = "U"
        gc.max_mdl_variants = [v]
        gc.max_mdl_reportable_variants = [v]
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "p.Thr34Ala"

    def test_s_shows_no_high_confidence(self, processor, make_lims_gene_code):
        gc = make_lims_gene_code()
        gc.max_mdl_interpretation = "S"
        gc.max_mdl_variants = []
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "No high confidence mutations detected"

    def test_wt_shows_no_mutations(self, processor, make_lims_gene_code):
        gc = make_lims_gene_code()
        gc.max_mdl_interpretation = "WT"
        gc.max_mdl_variants = []
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "No mutations detected"

    def test_na_shows_no_mutations(self, processor, make_lims_gene_code):
        gc = make_lims_gene_code()
        gc.max_mdl_interpretation = "NA"
        gc.max_mdl_variants = []
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "No mutations detected"

    def test_insufficient_coverage_shows_no_sequence(self, processor, make_lims_gene_code):
        gc = make_lims_gene_code()
        gc.max_mdl_interpretation = "Insufficient Coverage"
        gc.max_mdl_variants = []
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "No sequence"

    def test_rpob_synonymous_rrdr_appends_tag(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Phe433Phe", type="synonymous_variant",
        )
        gc.max_mdl_interpretation = "S"
        gc.max_mdl_variants = [v]
        gc.max_mdl_reportable_variants = [v]
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "p.Phe433Phe [synonymous]"

    def test_s_with_mixed_synonymous_only_reports_rrdr(self, processor, make_lims_gene_code, make_variant):
        """When max interpretation is S with both RRDR and non-RRDR synonymous variants, only the RRDR one is reported."""
        gc = make_lims_gene_code()
        v_syn_outside = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Ala400Ala", nucleotide_change="c.1200A>G",
            type="synonymous_variant",
        )
        v_syn_rrdr = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Phe433Phe", nucleotide_change="c.1299C>T",
            type="synonymous_variant",
        )
        gc.max_mdl_interpretation = "S"
        gc.max_mdl_variants = [v_syn_outside, v_syn_rrdr]
        gc.max_mdl_reportable_variants = [v_syn_rrdr]
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "p.Phe433Phe [synonymous]"
        assert "p.Ala400Ala" not in gc.gene_target_value

    def test_na_protein_change_falls_back_to_nt(self, processor, make_lims_gene_code, make_variant):
        """When protein_change is 'NA', gene_target should use nucleotide_change."""
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="R", nucleotide_change="c.1349C>T", protein_change="NA")
        gc.max_mdl_interpretation = "R"
        gc.max_mdl_variants = [v]
        gc.max_mdl_reportable_variants = [v]
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "c.1349C>T"


class TestResolveDrugTarget:
    def test_r_non_rpob_resistance(self, processor, make_lims_record, make_variant):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        v = make_variant(
            mdl_interpretation="R",
            gene_name="katG", drug="isoniazid", gene_id="Rv1908c", protein_change="p.Ser315Thr",
        )
        record.gene_codes["katG"].max_mdl_interpretation = "R"
        record.gene_codes["katG"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "Mutation(s) associated with resistance to isoniazid detected"

    def test_r_rpob_standard_mutation(self, processor, make_lims_record, make_variant):
        """Standard resistance R mutation is not affected by synonymous RRDR variant in reportable."""
        record = make_lims_record()
        v_r = make_variant(mdl_interpretation="R", gene_name="rpoB", drug="rifampicin", protein_change="p.Ser450Leu")
        v_syn_rrdr = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Phe433Phe", nucleotide_change="c.1299C>T",
            type="synonymous_variant",
        )
        record.gene_codes["rpoB"].max_mdl_interpretation = "R"
        record.gene_codes["rpoB"].max_mdl_variants = [v_r]
        record.gene_codes["rpoB"].max_mdl_reportable_variants = [v_r, v_syn_rrdr]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "Predicted resistance to rifampicin"

    def test_r_rpob_only_low_level_mutations(self, processor, make_lims_record, make_variant):
        """Low-level resistance R mutation is not affected by synonymous RRDR variant in reportable."""
        record = make_lims_record()
        v_r = make_variant(mdl_interpretation="R", gene_name="rpoB", drug="rifampicin", protein_change="p.His445Asn")
        v_syn_rrdr = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Phe433Phe", nucleotide_change="c.1299C>T",
            type="synonymous_variant",
        )
        record.gene_codes["rpoB"].max_mdl_interpretation = "R"
        record.gene_codes["rpoB"].max_mdl_variants = [v_r]
        record.gene_codes["rpoB"].max_mdl_reportable_variants = [v_r, v_syn_rrdr]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "Predicted low-level resistance to rifampicin. May test susceptible by phenotypic methods"

    def test_r_rpob_mixed_low_and_standard(self, processor, make_lims_record, make_variant):
        """When both low-level and standard rpoB mutations are present, standard wins."""
        record = make_lims_record()
        v1 = make_variant(
            mdl_interpretation="R",
            gene_name="rpoB", drug="rifampicin", protein_change="p.His445Asn", nucleotide_change="c.1333C>A",
        )
        v2 = make_variant(mdl_interpretation="R", gene_name="rpoB", drug="rifampicin", protein_change="p.Ser450Leu")
        record.gene_codes["rpoB"].max_mdl_interpretation = "R"
        record.gene_codes["rpoB"].max_mdl_variants = [v1, v2]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "Predicted resistance to rifampicin"

    def test_u_returns_uncertain_significance(self, processor, make_lims_record, make_variant):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        v = make_variant(mdl_interpretation="U", gene_name="katG", drug="isoniazid", gene_id="Rv1908c")
        record.gene_codes["katG"].max_mdl_interpretation = "U"
        record.gene_codes["katG"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "The detected mutation(s) have uncertain significance. Resistance to isoniazid cannot be ruled out"

    def test_s_returns_no_mutations(self, processor, make_lims_record, make_variant):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        v = make_variant(mdl_interpretation="S", gene_name="katG", drug="isoniazid", gene_id="Rv1908c")
        record.gene_codes["katG"].max_mdl_interpretation = "S"
        record.gene_codes["katG"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "No mutations associated with resistance to isoniazid detected"

    def test_na_returns_no_mutations(self, processor, make_lims_record):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        record.gene_codes["katG"].max_mdl_interpretation = "NA"
        record.gene_codes["katG"].max_mdl_variants = []
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "No mutations associated with resistance to isoniazid detected"

    def test_rpob_s_with_synonymous_rrdr(self, processor, make_lims_record, make_variant):
        record = make_lims_record()
        v = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampicin",
            protein_change="p.Phe433Phe", type="synonymous_variant",
        )
        record.gene_codes["rpoB"].max_mdl_interpretation = "S"
        record.gene_codes["rpoB"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "Predicted susceptibility to rifampicin. The detected synonymous mutation(s) do not confer resistance"

    def test_rpob_s_without_synonymous_rrdr(self, processor, make_lims_record, make_variant):
        record = make_lims_record()
        v = make_variant(mdl_interpretation="S", gene_name="rpoB", drug="rifampicin", protein_change="p.Ser450Leu")
        record.gene_codes["rpoB"].max_mdl_interpretation = "S"
        record.gene_codes["rpoB"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert record.drug_target_value == "Predicted susceptibility to rifampicin"


class TestProcessLimsMtbcId:
    def test_get_pnca_his57asp_variants(self, processor, make_variant):
        v1 = make_variant(gene_name="pncA", drug="pyrazinamide", protein_change="p.His57Asp")
        v2 = make_variant(gene_name="pncA", drug="pyrazinamide", protein_change="p.Ser104Arg")
        v3 = make_variant(gene_name="rpoB", drug="rifampicin", protein_change="p.Ser450Leu")
        variants = [v1, v2, v3]
        pnca_variants = processor._get_pnca_his57asp_variants(variants)
        assert pnca_variants == [v1]

    def test_wgs_lineage4_detected(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "lineage4", "lineage4.1"
        )
        assert "DNA of Mycobacterium tuberculosis species detected" == result

    def test_wgs_bcg_detected(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "BCG", "BCG"
        )
        assert "DNA of Mycobacterium bovis BCG detected" == result

    def test_wgs_bovis_la1_detected(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "La1", "La1.1"
        )
        assert "DNA of Mycobacterium bovis (not BCG) detected" == result

    def test_wgs_empty_lineage_falls_back_to_complex(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "", ""
        )
        assert "DNA of Mycobacterium tuberculosis complex detected" == result

    def test_wgs_na_lineage_falls_back_to_complex(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "NA", ""
        )
        assert "DNA of Mycobacterium tuberculosis complex detected" == result

    def test_tngs_pnca_his57asp_passes_qc(self, processor, make_variant, make_lims_record, make_lims_gene_code, make_locus_coverage, mock_config):
        lims_records = [
            make_lims_record(
                drug="pyrazinamide",
                drug_code="PZA",
                gene_codes={"pncA": make_lims_gene_code("PZA_pncA")}
            )
        ]
        locus_coverage_map = {"Rv2043c": make_locus_coverage(locus_tag="Rv2043c", breadth_of_coverage=1.0)}

        mock_config.TNGS = True
        v = make_variant(
            gene_name="pncA", drug="pyrazinamide",
            protein_change="p.His57Asp", gene_id="Rv2043c",
        )
        v.fails_positional_qc = False
        result = processor.process_lims_mtbc_id(
            lims_records, [v], locus_coverage_map, "", ""
        )
        assert "DNA of Mycobacterium bovis detected" == result

    def test_tngs_pnca_his57asp_ignores_fails_qc(self, processor, make_variant, make_lims_record, make_lims_gene_code, make_locus_coverage, mock_config):
        lims_records = [
            make_lims_record(
                drug="pyrazinamide",
                drug_code="PZA",
                gene_codes={"pncA" : make_lims_gene_code("PZA_pncA")}
            )
        ]
        locus_coverage_map = {"Rv2043c": make_locus_coverage(locus_tag="Rv2043c", breadth_of_coverage=1.0)}

        mock_config.TNGS = True
        v = make_variant(
            gene_name="pncA", drug="pyrazinamide",
            protein_change="p.His57Asp", gene_id="Rv2043c",
        )
        v.fails_positional_qc = True
        result = processor.process_lims_mtbc_id(
            lims_records, [v], locus_coverage_map, "", ""
        )
        assert "DNA of Mycobacterium tuberculosis complex detected (M. bovis not ruled out)" == result

    def test_tngs_no_pnca_his57asp(self, processor, make_variant, make_lims_record, make_lims_gene_code, make_locus_coverage, mock_config):
        lims_records = [
            make_lims_record(
                drug="pyrazinamide",
                drug_code="PZA",
                gene_codes={"pncA" : make_lims_gene_code("PZA_pncA")}
            )
        ]
        locus_coverage_map = {"Rv2043c": make_locus_coverage(locus_tag="Rv2043c", breadth_of_coverage=1.0)}

        mock_config.TNGS = True
        v = make_variant(
            gene_name="pncA", drug="pyrazinamide",
            protein_change="p.His59Asp", gene_id="Rv2043c",
        )
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "", ""
        )
        assert "DNA of Mycobacterium tuberculosis complex detected (not M. bovis)" == result

    def test_low_coverage_returns_not_detected(self, processor, make_lims_record, make_locus_coverage, mock_config):
        mock_config.MIN_PERCENT_LOCI_COVERED = 1.0  # require 100%
        lims_records = [make_lims_record()]
        low_cov_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], low_cov_map, "lineage4", "lineage4.1"
        )
        assert "DNA of Mycobacterium tuberculosis complex NOT detected" == result


class TestPassesLimsCoverageFractionERR:
    """Tests for _passes_lims_coverage_fraction with USE_ERR_AS_BRR flag."""

    def _make_err(self, breadth=0.95, coords=None):
        if coords is None:
            coords = [(100, 150), (250, 350)]
        return ERRCoverage(coords=coords, breadth_of_coverage=breadth, average_depth=100.0)

    def test_err_high_boc_passes_lims_qc(self, processor, make_lims_record, make_locus_coverage, mock_config):
        mock_config.USE_ERR_AS_BRR = True
        lims_records = [make_lims_record()]
        locus_coverage_map = {
            "Rv0667": make_locus_coverage(
                locus_tag="Rv0667", breadth_of_coverage=0.50,
                err_coverage=self._make_err(breadth=0.95),
            )
        }
        assert processor._passes_lims_coverage_fraction(lims_records, locus_coverage_map) is True

    def test_err_low_boc_fails_lims_qc(self, processor, make_lims_record, make_locus_coverage, mock_config):
        mock_config.USE_ERR_AS_BRR = True
        mock_config.MIN_PERCENT_LOCI_COVERED = 1.0  # require 100%
        lims_records = [make_lims_record()]
        err_coverage = self._make_err(breadth=0.50)
        locus_coverage_map = {
            "Rv0667": make_locus_coverage(
                locus_tag="Rv0667", breadth_of_coverage=0.50,
                err_coverage=err_coverage,
            )
        }
        assert processor._passes_lims_coverage_fraction(lims_records, locus_coverage_map) is False

    def test_err_valid_deletion_overrides_low_boc(self, processor, make_lims_record, make_locus_coverage, make_variant, mock_config):
        mock_config.USE_ERR_AS_BRR = True
        mock_config.MIN_PERCENT_LOCI_COVERED = 1.0
        del_variant = make_variant(nucleotide_change="c.1_100del", gene_id="Rv0667", pos=120)

        err = self._make_err(breadth=0.50, coords=[(100, 150), (250, 350)])
        err.valid_deletions.append(del_variant)

        lims_records = [make_lims_record()]
        locus_coverage_map = {
            "Rv0667": make_locus_coverage(
                locus_tag="Rv0667", breadth_of_coverage=0.50,
                err_coverage=err,
            )
        }
        assert processor._passes_lims_coverage_fraction(lims_records, locus_coverage_map) is True

    def test_flag_off_falls_back_to_locus(self, processor, make_lims_record, make_locus_coverage, mock_config):
        mock_config.USE_ERR_AS_BRR = False
        mock_config.MIN_PERCENT_LOCI_COVERED = 1.0
        lims_records = [make_lims_record()]
        locus_coverage_map = {
            "Rv0667": make_locus_coverage(
                locus_tag="Rv0667", breadth_of_coverage=0.50,
                err_coverage=self._make_err(breadth=0.95),
            )
        }
        # Flag off -> uses locus breadth (0.50) -> fails
        assert processor._passes_lims_coverage_fraction(lims_records, locus_coverage_map) is False
