import logging

import yaml

from pathlib import Path

from tbp_parser.LIMS.lims_yml_builder import (
    DRUG_CODES,
    build_lims_fmt,
    build_lims_report_format,
    write_lims_report_format_yml,
)
from tbp_parser.LIMS.lims_yml_parser import parse_lims_yml_file

TEST_GENE_DB = str(Path(__file__).parent / "test_files" / "test_gene_db.yml")


def _write_gene_database(tmp_path, gene_database):
    """Write a gene database dict to a YAML file and return its path, since the builder reads from disk."""
    gene_database_yml = tmp_path / "gene-database.yml"
    gene_database_yml.write_text(yaml.dump(gene_database))
    return str(gene_database_yml)


class TestBuildLimsReportFormat:

    def test_inverts_a_gene_into_every_drug_it_is_associated_with(self, tmp_path):
        # the gene database is keyed by gene, the LIMS format by drug, so a gene with several
        # drugs has to appear once under each of them, with a per-drug gene code.
        # the full-equality assertion is also this file's baseline for the shape of an entry.
        gene_database_yml = _write_gene_database(tmp_path, {
            "Rv0006": {"gene_name": "gyrA", "drugs": ["levofloxacin", "moxifloxacin"]},
        })

        lims_report_format = build_lims_report_format(gene_database_yml)

        assert lims_report_format == [
            {"drug": "levofloxacin", "drug_code": "LFX", "gene_codes": {"gyrA": "LFX_gyrA"}},
            {"drug": "moxifloxacin", "drug_code": "MFX", "gene_codes": {"gyrA": "MFX_gyrA"}},
        ]

    def test_sorts_drugs_alphabetically_and_genes_case_insensitively(self, tmp_path):
        # deliberately unsorted, and mixed case, so this fails if either sort is dropped.
        # a case-sensitive sort would put the capitalized Rv2477c ahead of every lowercase gene.
        gene_database_yml = _write_gene_database(tmp_path, {
            "Rv2477c": {"gene_name": "Rv2477c", "drugs": ["kanamycin", "amikacin"]},
            "Rv1694": {"gene_name": "tlyA", "drugs": ["amikacin"]},
            "EBG00000313325": {"gene_name": "rrs", "drugs": ["amikacin"]},
        })

        lims_report_format = build_lims_report_format(gene_database_yml)

        assert [entry["drug"] for entry in lims_report_format] == ["amikacin", "kanamycin"]
        assert list(lims_report_format[0]["gene_codes"]) == ["rrs", "Rv2477c", "tlyA"]

    def test_skips_genes_without_drugs(self, tmp_path, caplog):
        # a gene associated with no drug has no drug column to sit under, so it cannot be reported
        gene_database_yml = _write_gene_database(tmp_path, {
            "Rv0667": {"gene_name": "rpoB", "drugs": ["rifampicin"]},
            "Rv0001": {"gene_name": "dnaA", "drugs": []},
            "Rv0002": {"gene_name": "dnaN"},
        })

        with caplog.at_level(logging.WARNING, logger="tbp_parser.LIMS.lims_yml_builder"):
            lims_report_format = build_lims_report_format(gene_database_yml)

        assert [entry["drug"] for entry in lims_report_format] == ["rifampicin"]
        assert any(
            r.levelno == logging.WARNING and "2 genes have no associated drugs" in r.message and "dnaA, dnaN" in r.message
            for r in caplog.records
        )

    def test_falls_back_to_the_locus_tag_when_an_entry_has_no_gene_name(self, tmp_path):
        gene_database_yml = _write_gene_database(tmp_path, {
            "Rv0010c": {"drugs": ["isoniazid"]},
        })

        lims_report_format = build_lims_report_format(gene_database_yml)

        assert lims_report_format[0]["gene_codes"] == {"Rv0010c": "INH_Rv0010c"}

    def test_derives_a_drug_code_for_an_unknown_drug(self, tmp_path, caplog):
        # an unrecognized drug still gets a column rather than being dropped, but it is warned
        # about because the derived code is a guess that LIMS may not accept
        gene_database_yml = _write_gene_database(tmp_path, {
            "Rv0667": {"gene_name": "rpoB", "drugs": ["novelcin"]},
        })

        with caplog.at_level(logging.WARNING, logger="tbp_parser.LIMS.lims_yml_builder"):
            lims_report_format = build_lims_report_format(gene_database_yml)

        assert lims_report_format == [
            {"drug": "novelcin", "drug_code": "NOVELCIN", "gene_codes": {"rpoB": "NOVELCIN_rpoB"}}
        ]
        assert any(
            r.levelno == logging.WARNING and "novelcin (NOVELCIN)" in r.message
            for r in caplog.records
        )


class TestWriteLimsReportFormatYml:

    def test_round_trips_through_the_lims_report_format(self, tmp_path):
        output_path = tmp_path / "lims-report-format.yml"

        lims_report_format = build_lims_report_format(TEST_GENE_DB)
        write_lims_report_format_yml(lims_report_format, str(output_path))

        assert yaml.safe_load(output_path.read_text()) == lims_report_format

    def test_preserves_key_order_within_an_entry(self, tmp_path):
        # sort_keys=False keeps each entry readable in the written file, drug first
        output_path = tmp_path / "lims-report-format.yml"

        lims_report_format = build_lims_report_format(TEST_GENE_DB)
        write_lims_report_format_yml(lims_report_format, str(output_path))

        assert list(yaml.safe_load(output_path.read_text())[0]) == ["drug", "drug_code", "gene_codes"]


class TestBuildLimsFmt:

    def test_writes_a_format_the_parser_can_read_back(self, tmp_path):
        # the whole point of the subcommand is producing a --lims_report_format_yml file, so the
        # written file has to load as LIMSRecords. this is also the only pass over the real gene
        # database, so it checks that no gene/drug pair is dropped along the way.
        gene_database = yaml.safe_load(Path(TEST_GENE_DB).read_text())
        expected_gene_code_count = sum(len(entry["drugs"]) for entry in gene_database.values())
        output_path = tmp_path / "lims-report-format.yml"

        build_lims_fmt(TEST_GENE_DB, str(output_path))
        lims_records = parse_lims_yml_file(str(output_path))

        assert [record.drug for record in lims_records] == sorted(DRUG_CODES)
        assert sum(len(record.gene_codes) for record in lims_records) == expected_gene_code_count
        rifampicin = next(record for record in lims_records if record.drug == "rifampicin")
        assert rifampicin.drug_code == "RIF"
        assert rifampicin.gene_codes["rpoB"].gene_code == "RIF_rpoB"
