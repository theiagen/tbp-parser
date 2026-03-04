import pandas as pd
from pathlib import Path

from Reporters.lab_report import write_laboratorian_report
from Reporters.lims_report import write_lims_report
from Reporters.looker_report import write_looker_report
from Reporters.coverage_report import write_coverage_report
from Coverage.coverage_data import ERRCoverage
from LIMS import LIMSRecord, LIMSGeneCode


def _run_report(mock_config, tmp_path, writer_fn, report_suffix, *args):
    """Set output prefix, run a report writer, assert the file exists, and return the DataFrame."""
    mock_config.OUTPUT_PREFIX = str(tmp_path / "test")
    writer_fn(*args)
    output = Path(f"{mock_config.OUTPUT_PREFIX}.{report_suffix}.csv")
    assert output.exists()
    return pd.read_csv(output)


class TestWriteLaboratorianReport:
    def test_writes_csv(self, mock_config, make_variant, tmp_path):
        v = make_variant(mdl_interpretation="R", rationale="WHO classification")
        df = _run_report(mock_config, tmp_path, write_laboratorian_report, "laboratorian_report", [v])
        assert len(df) == 1
        assert df["sample_id"].iloc[0] == "test_sample"
        assert df["tbprofiler_gene_name"].iloc[0] == "rpoB"
        # Variant normalizes "rifampicin" -> "rifampin"
        assert df["antimicrobial"].iloc[0] == "rifampin"
        assert df["mdl_interpretation"].iloc[0] == "R"

    def test_multiple_variants(self, mock_config, make_variant, tmp_path):
        v1 = make_variant(gene_name="rpoB", drug="rifampicin", mdl_interpretation="R", rationale="WHO classification")
        v2 = make_variant(gene_name="katG", gene_id="Rv1908c", drug="isoniazid", mdl_interpretation="R", rationale="WHO classification")
        df = _run_report(mock_config, tmp_path, write_laboratorian_report, "laboratorian_report", [v1, v2])
        assert len(df) == 2

    def test_warning_joined(self, mock_config, make_variant, tmp_path):
        v = make_variant(mdl_interpretation="Insufficient Coverage", rationale="NA")
        v.warning = {"Failed quality in the mutation position", "Insufficient coverage in locus"}
        df = _run_report(mock_config, tmp_path, write_laboratorian_report, "laboratorian_report", [v])
        warning_val = df["warning"].iloc[0]
        assert "Failed quality" in warning_val
        assert "Insufficient coverage" in warning_val


class TestWriteLimsReport:
    def test_writes_lims_and_transposed(self, mock_config, tmp_path):
        gc = LIMSGeneCode(gene_code="M_DST_D02_rpoB")
        gc.gene_target_value = "p.Ser450Leu"
        rec = LIMSRecord(
            drug="rifampin", drug_code="M_DST_D02",
            gene_codes={"rpoB": gc},
        )
        rec.drug_target_value = "Predicted resistance to rifampin"

        mock_config.OUTPUT_PREFIX = str(tmp_path / "test")
        write_lims_report([rec], "DNA of M. tuberculosis", "test_sample", "lineage4")

        output = Path(f"{mock_config.OUTPUT_PREFIX}.lims_report.csv")
        transposed = Path(f"{mock_config.OUTPUT_PREFIX}.lims_report.transposed.csv")
        assert output.exists()
        assert transposed.exists()

        df = pd.read_csv(output)
        assert "M_DST_D02" in df.columns
        assert "Predicted resistance" in df["M_DST_D02"].iloc[0]

    def test_gene_target_value_written(self, mock_config, tmp_path):
        gc = LIMSGeneCode(gene_code="INH_katG")
        gc.gene_target_value = "p.Ser315Thr"
        rec = LIMSRecord(
            drug="isoniazid", drug_code="INH",
            gene_codes={"katG": gc},
        )
        rec.drug_target_value = "Mutation(s) associated with resistance to isoniazid detected"

        df = _run_report(mock_config, tmp_path, write_lims_report, "lims_report", [rec], "DNA of M. tuberculosis", "test_sample", "lineage4")
        assert "p.Ser315Thr" in df["INH_katG"].iloc[0]

    def test_lineage_written_to_metadata_columns(self, mock_config, tmp_path):
        gc = LIMSGeneCode(gene_code="M_DST_D02_rpoB")
        gc.gene_target_value = "No mutations detected"
        rec = LIMSRecord(
            drug="rifampin", drug_code="M_DST_D02",
            gene_codes={"rpoB": gc},
        )
        rec.drug_target_value = "No mutations detected"

        df = _run_report(mock_config, tmp_path, write_lims_report, "lims_report", [rec], "DNA of M. tuberculosis", "test_sample", "lineage4")
        assert df["M_DST_A01_ID"].iloc[0] == "DNA of M. tuberculosis"
        assert df["M_DST_O01_Lineage"].iloc[0] == "lineage4"
        assert df["MDL sample accession numbers"].iloc[0] == "test_sample"


class TestWriteLookerReport:
    def test_writes_looker_csv(self, mock_config, make_variant, tmp_path):
        v = make_variant(gene_name="rpoB", drug="rifampicin")
        v.looker_interpretation = "R"
        df = _run_report(mock_config, tmp_path, write_looker_report, "looker_report", [v], "DNA of Mycobacterium tuberculosis species detected", "test_sample", "lineage4")
        # Variant normalizes "rifampicin" -> "rifampin"
        assert df["rifampin"].iloc[0] == "R"
        assert df["lineage"].iloc[0] == "lineage4"

    def test_highest_interpretation_wins_for_same_drug(self, mock_config, make_variant, tmp_path):
        v1 = make_variant(gene_name="rpoB", drug="rifampicin", nucleotide_change="c.1349C>T",
                          protein_change="p.Ser450Leu")
        v1.looker_interpretation = "R"
        v2 = make_variant(gene_name="rpoB", drug="rifampicin", nucleotide_change="c.1300A>G",
                          protein_change="p.Thr434Ala")
        v2.looker_interpretation = "S"
        df = _run_report(mock_config, tmp_path, write_looker_report, "looker_report", [v1, v2], "lineage", "test_sample", "lineage4")
        assert df["rifampin"].iloc[0] == "R"

    def test_sample_id_and_operator_written(self, mock_config, make_variant, tmp_path):
        mock_config.OPERATOR = "my_operator"
        v = make_variant(gene_name="rpoB", drug="rifampicin")
        v.looker_interpretation = "S"
        df = _run_report(mock_config, tmp_path, write_looker_report, "looker_report", [v], "lineage_text", "test_sample", "lineage4")
        assert df["sample_id"].iloc[0] == "test_sample"
        assert df["operator"].iloc[0] == "my_operator"
        assert df["ID"].iloc[0] == "lineage_text"


class TestWriteCoverageReports:
    def test_target_coverage_report(self, mock_config, make_target_coverage, tmp_path):
        gene_map = {
            "rpoB": make_target_coverage(gene_name="rpoB", breadth_of_coverage=0.995, average_depth=150.0),
            "katG": make_target_coverage(gene_name="katG", locus_tag="Rv1908c", breadth_of_coverage=0.85, average_depth=50.0),
        }
        df = _run_report(mock_config, tmp_path, write_coverage_report, "target_coverage_report", "test_sample", gene_map)
        assert len(df) == 2
        assert "gene_name" in df.columns
        assert "percent_coverage" in df.columns
        assert "average_depth" in df.columns

    def test_locus_coverage_report(self, mock_config, make_locus_coverage, tmp_path):
        locus_map = {
            "Rv0667": make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.995, average_depth=150.0),
        }
        df = _run_report(mock_config, tmp_path, write_coverage_report, "locus_coverage_report", "test_sample", locus_map)
        assert len(df) == 1
        assert "locus_tag" in df.columns

    def test_target_report_uses_err_deletions(self, mock_config, make_target_coverage, make_variant, tmp_path):
        del_variant = make_variant(nucleotide_change="c.1_100del", pos=120)
        err = ERRCoverage(coords=[(100, 150)], breadth_of_coverage=0.95, average_depth=100.0, valid_deletions=[del_variant])
        tc = make_target_coverage(gene_name="rpoB", err_coverage=err)

        df = _run_report(mock_config, tmp_path, write_coverage_report, "target_coverage_report", "test_sample", {"rpoB": tc})
        assert "c.1_100del" in df["qc_warning"].iloc[0]

    def test_target_report_no_err_uses_target_deletions(self, mock_config, make_target_coverage, make_variant, tmp_path):
        del_variant = make_variant(nucleotide_change="c.200_300del", pos=120)
        tc = make_target_coverage(gene_name="rpoB", err_coverage=None)
        tc.valid_deletions.append(del_variant)

        df = _run_report(mock_config, tmp_path, write_coverage_report, "target_coverage_report", "test_sample", {"rpoB": tc})
        assert "c.200_300del" in df["qc_warning"].iloc[0]

    def test_locus_report_uses_err_deletions(self, mock_config, make_locus_coverage, make_variant, tmp_path):
        del_variant = make_variant(nucleotide_change="c.1_100del", pos=120)
        err = ERRCoverage(coords=[(100, 150), (250, 350)], breadth_of_coverage=0.95, average_depth=100.0, valid_deletions=[del_variant])
        lc = make_locus_coverage(locus_tag="Rv0667", err_coverage=err)

        df = _run_report(mock_config, tmp_path, write_coverage_report, "locus_coverage_report", "test_sample", {"Rv0667": lc})
        assert "c.1_100del" in df["qc_warning"].iloc[0]
