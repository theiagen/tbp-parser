import pytest
from LIMS import LIMSProcessor, LIMSRecord, LIMSGeneCode
from Coverage.coverage_data import ERRCoverage

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
    def _make(drug="rifampin", drug_code="M_DST_D02", gene_codes=None):
        if gene_codes is None:
            gene_codes = {"rpoB": LIMSGeneCode(gene_code="M_DST_D02_rpoB")}
        return LIMSRecord(drug=drug, drug_code=drug_code, gene_codes=gene_codes)
    return _make

class TestResistanceRanking:
    @pytest.mark.parametrize("higher,lower", [
        ("R", "U"),
        ("U", "S"),
        ("S", "WT"),
        ("WT", "Insufficient Coverage"),
        ("Insufficient Coverage", "NA"),
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


class TestResolveGeneTarget:
    def test_r_shows_protein_change(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="R", protein_change="p.Ser450Leu")
        gc.max_mdl_interpretation = "R"
        gc.max_mdl_variants = [v]
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "p.Ser450Leu"

    def test_u_shows_protein_change(self, processor, make_lims_gene_code, make_variant):
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="U", protein_change="p.Thr34Ala")
        gc.max_mdl_interpretation = "U"
        gc.max_mdl_variants = [v]
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

    def test_na_shows_no_sequence(self, processor, make_lims_gene_code):
        gc = make_lims_gene_code()
        gc.max_mdl_interpretation = "NA"
        gc.max_mdl_variants = []
        processor.resolve_gene_target(gc)
        assert gc.gene_target_value == "No sequence"

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
            gene_name="rpoB", drug="rifampin",
            protein_change="p.Phe433Phe", type="synonymous_variant",
        )
        gc.max_mdl_interpretation = "S"
        gc.max_mdl_variants = [v]
        processor.resolve_gene_target(gc)
        assert "[synonymous]" in gc.gene_target_value
        assert "p.Phe433Phe" in gc.gene_target_value

    def test_na_protein_change_falls_back_to_nt(self, processor, make_lims_gene_code, make_variant):
        """When protein_change is 'NA', gene_target should use nucleotide_change."""
        gc = make_lims_gene_code()
        v = make_variant(mdl_interpretation="R", nucleotide_change="c.1349C>T", protein_change="NA")
        gc.max_mdl_interpretation = "R"
        gc.max_mdl_variants = [v]
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
        assert "resistance to isoniazid detected" in record.drug_target_value

    def test_r_rpob_standard_mutation(self, processor, make_lims_record, make_variant):
        record = make_lims_record()
        v = make_variant(mdl_interpretation="R", gene_name="rpoB", drug="rifampin", protein_change="p.Ser450Leu")
        record.gene_codes["rpoB"].max_mdl_interpretation = "R"
        record.gene_codes["rpoB"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert "Predicted resistance to rifampin" in record.drug_target_value

    def test_r_rpob_only_low_level_mutations(self, processor, make_lims_record, make_variant):
        record = make_lims_record()
        v = make_variant(mdl_interpretation="R", gene_name="rpoB", drug="rifampin", protein_change="p.His445Asn")
        record.gene_codes["rpoB"].max_mdl_interpretation = "R"
        record.gene_codes["rpoB"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert "low-level resistance" in record.drug_target_value

    def test_r_rpob_mixed_low_and_standard(self, processor, make_lims_record, make_variant):
        """When both low-level and standard rpoB mutations are present, standard wins."""
        record = make_lims_record()
        v1 = make_variant(
            mdl_interpretation="R",
            gene_name="rpoB", drug="rifampin", protein_change="p.His445Asn", nucleotide_change="c.1333C>A",
        )
        v2 = make_variant(mdl_interpretation="R", gene_name="rpoB", drug="rifampin", protein_change="p.Ser450Leu")
        record.gene_codes["rpoB"].max_mdl_interpretation = "R"
        record.gene_codes["rpoB"].max_mdl_variants = [v1, v2]
        processor.resolve_drug_target(record)
        assert "Predicted resistance to rifampin" in record.drug_target_value

    def test_u_returns_uncertain_significance(self, processor, make_lims_record, make_variant):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        v = make_variant(mdl_interpretation="U", gene_name="katG", drug="isoniazid", gene_id="Rv1908c")
        record.gene_codes["katG"].max_mdl_interpretation = "U"
        record.gene_codes["katG"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert "uncertain significance" in record.drug_target_value

    def test_s_returns_no_mutations(self, processor, make_lims_record, make_variant):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        v = make_variant(mdl_interpretation="S", gene_name="katG", drug="isoniazid", gene_id="Rv1908c")
        record.gene_codes["katG"].max_mdl_interpretation = "S"
        record.gene_codes["katG"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert "No mutations" in record.drug_target_value

    def test_na_returns_pending_retest(self, processor, make_lims_record):
        record = make_lims_record(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": LIMSGeneCode(gene_code="INH_katG")},
        )
        record.gene_codes["katG"].max_mdl_interpretation = "NA"
        record.gene_codes["katG"].max_mdl_variants = []
        processor.resolve_drug_target(record)
        assert "Pending Retest" in record.drug_target_value

    def test_rpob_s_with_synonymous_rrdr(self, processor, make_lims_record, make_variant):
        record = make_lims_record()
        v = make_variant(
            mdl_interpretation="S",
            gene_name="rpoB", drug="rifampin",
            protein_change="p.Phe433Phe", type="synonymous_variant",
        )
        record.gene_codes["rpoB"].max_mdl_interpretation = "S"
        record.gene_codes["rpoB"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert "synonymous mutation" in record.drug_target_value

    def test_rpob_s_without_synonymous_rrdr(self, processor, make_lims_record, make_variant):
        record = make_lims_record()
        v = make_variant(mdl_interpretation="S", gene_name="rpoB", drug="rifampin", protein_change="p.Ser450Leu")
        record.gene_codes["rpoB"].max_mdl_interpretation = "S"
        record.gene_codes["rpoB"].max_mdl_variants = [v]
        processor.resolve_drug_target(record)
        assert "Predicted susceptibility to rifampin" in record.drug_target_value
        assert "synonymous" not in record.drug_target_value


class TestProcessLimsMtbcId:
    def test_wgs_lineage4_detected(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "lineage4", "lineage4.1"
        )
        assert "Mycobacterium tuberculosis species detected" in result

    def test_wgs_bcg_detected(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "BCG", "BCG"
        )
        assert "bovis BCG" in result

    def test_wgs_bovis_la1_detected(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "La1", "La1.1"
        )
        assert "bovis" in result
        assert "not BCG" in result

    def test_wgs_empty_lineage_falls_back_to_complex(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "", ""
        )
        assert "tuberculosis complex detected" in result
        assert "NOT" not in result

    def test_wgs_na_lineage_falls_back_to_complex(self, processor, make_lims_record, make_locus_coverage):
        lims_records = [make_lims_record()]
        locus_coverage_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=1.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], locus_coverage_map, "NA", ""
        )
        assert "tuberculosis complex detected" in result

    def test_tngs_pnca_his57asp_passes_qc(
        self, processor, make_variant, make_lims_record, make_lims_gene_code, make_locus_coverage, mock_config
    ):
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
        v.warning = set()
        result = processor.process_lims_mtbc_id(
            lims_records, [v], locus_coverage_map, "", ""
        )
        assert "Mycobacterium bovis detected" in result
        assert "not ruled out" not in result

    def test_tngs_pnca_his57asp_ignores_fails_qc(
        self, processor, make_variant, make_lims_record, make_lims_gene_code, make_locus_coverage, mock_config
    ):
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
        v.warning = {"Failed quality in the mutation position"}
        result = processor.process_lims_mtbc_id(
            lims_records, [v], locus_coverage_map, "", ""
        )
        assert "M. bovis not ruled out" in result

    def test_tngs_no_pnca_his57asp(
        self, processor, make_variant, make_lims_record, make_lims_gene_code, make_locus_coverage, mock_config
    ):
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
        assert "not M. bovis" in result

    def test_low_coverage_returns_not_detected(
        self, processor, make_lims_record, make_locus_coverage, mock_config
    ):
        mock_config.MIN_PERCENT_LOCI_COVERED = 1.0  # require 100%
        lims_records = [make_lims_record()]
        low_cov_map = {"Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.0)}
        result = processor.process_lims_mtbc_id(
            lims_records, [], low_cov_map, "lineage4", "lineage4.1"
        )
        assert "NOT detected" in result


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
