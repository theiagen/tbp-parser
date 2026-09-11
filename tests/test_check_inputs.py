import logging

import pytest

from tbp_parser.Utilities.check_inputs import validate_inputs


class TestValidateInputs:
    """
    Tests for validating the BED, LIMS, and results JSON input files against the gene database.

    Every test here runs against the gene database that `setup_gene_database` in conftest.py loads
    from `tests/test_files/test_gene_db.yml`. That fixture is session-scoped and autouse=True, so it is
    already in place without being requested.
    """

    def test_logs_success_when_all_inputs_are_consistent(self, make_bed_record, make_lims_record, caplog):
        # Asserting on the success log, rather than just on nothing being raised.
        # Separates every check ran and passed" from "no check ran at all". An empty variant_records list is
        # the no-results-JSON case, which is skipped rather than failed.
        bed_records = [make_bed_record(locus_tag="Rv0006", gene_name="gyrA")]
        lims_records = [make_lims_record("levofloxacin", "LFX", "gyrA")]

        with caplog.at_level(logging.INFO, logger="tbp_parser.Utilities.check_inputs"):
            validate_inputs(bed_records, lims_records, [])

        assert any(
            r.levelno == logging.INFO and "are present in the Gene Database" in r.message
            for r in caplog.records
        )

    def test_accepts_a_subset_of_the_gene_database(self, make_bed_record):
        # fewer genes and drugs than the gene database contains is expected and allowed
        validate_inputs([make_bed_record(locus_tag="Rv0006", gene_name="gyrA")], [], [])

    def test_rejects_bed_locus_tag_missing_from_database(self, make_bed_record):
        bed_records = [make_bed_record(locus_tag="RvNotAGene", gene_name="nope")]

        with pytest.raises(ValueError, match="genes from the BED file are missing in the Gene Database: RvNotAGene"):
            validate_inputs(bed_records, [], [])

    def test_rejects_lims_gene_drug_pair_missing_from_database(self, make_bed_record, make_lims_record):
        # gyrA is a real gene, but it is not associated with bedaquiline
        lims_records = [make_lims_record("bedaquiline", "BDQ", "gyrA")]

        with pytest.raises(ValueError, match=r"missing in the Gene Database: bedaquiline\|gyrA\|Rv0006"):
            validate_inputs([make_bed_record(locus_tag="Rv0006", gene_name="gyrA")], lims_records, [])

    def test_rejects_lims_drug_missing_from_database(self, make_bed_record, make_lims_record):
        # a drug no gene in the gene database is associated with is reported once, on its own, rather
        # than once per gene listed under it
        lims_records = [make_lims_record("notADrug", "NAD", "gyrA", "gyrB")]

        with pytest.raises(ValueError) as excinfo:
            validate_inputs(
                [
                    make_bed_record(locus_tag="Rv0006", gene_name="gyrA"),
                    make_bed_record(locus_tag="Rv0005", gene_name="gyrB"),
                ],
                lims_records,
                [],
            )

        message = str(excinfo.value)
        assert "drugs from the LIMS report format yaml file are missing in the Gene Database: notADrug" in message
        assert "notADrug|gyrA|Rv0006" not in message
        assert "notADrug|gyrB|Rv0005" not in message

    def test_rejects_lims_gene_missing_from_the_bed_file(self, make_bed_record, make_lims_record):
        lims_records = [make_lims_record("levofloxacin", "LFX", "gyrA")]

        with pytest.raises(ValueError, match=r"missing in the BED file: gyrA\|Rv0006"):
            validate_inputs([make_bed_record(locus_tag="Rv0005", gene_name="gyrB")], lims_records, [])

    def test_resolves_lims_genes_by_locus_tag_and_alias(self, make_bed_record, make_lims_record):
        # every getter resolves gene names, locus tags, and aliases interchangeably
        bed_records = [
            make_bed_record(locus_tag="Rv0006", gene_name="gyrA"),
            make_bed_record(locus_tag="EBG00000313325", gene_name="rrs"),
        ]
        lims_records = [
            make_lims_record("levofloxacin", "LFX", "Rv0006"),
            make_lims_record("amikacin", "AMK", "Rvnr01"),
        ]

        validate_inputs(bed_records, lims_records, [])

    def test_reports_every_failure_in_a_single_error(self, make_bed_record, make_lims_record):
        bed_records = [make_bed_record(locus_tag="RvNotAGene", gene_name="nope")]
        lims_records = [
            make_lims_record("levofloxacin", "LFX", "alsoNotAGene"),
            make_lims_record("bedaquiline", "BDQ", "gyrA"),
        ]

        with pytest.raises(ValueError) as excinfo:
            validate_inputs(bed_records, lims_records, [])

        message = str(excinfo.value)
        assert "RvNotAGene" in message
        assert "alsoNotAGene" in message
        assert "bedaquiline|gyrA|Rv0006" in message
        assert "build_gene_db" in message

    def test_reports_an_unresolvable_lims_gene_as_unchecked_against_the_bed_file(self, make_bed_record, make_lims_record):
        # with no locus tag there is nothing to compare, and the BED file may well carry a region for
        # the gene under a locus tag the gene database does not know, so it must not be called missing
        lims_records = [make_lims_record("levofloxacin", "LFX", "notAGene")]

        with pytest.raises(ValueError) as excinfo:
            validate_inputs([make_bed_record(locus_tag="Rv0006", gene_name="gyrA")], lims_records, [])

        message = str(excinfo.value)
        assert "missing in the Gene Database: notAGene" in message
        assert "could not be checked against the BED file" in message
        assert "missing in the BED file" not in message


class TestValidateInputsAgainstResultsJson:
    """Tests for validating the results JSON gene/drug pairings against the gene database."""

    def test_accepts_pairings_present_in_the_gene_database(self, make_bed_record, make_variant_record, make_annotation):
        variant_records = [
            make_variant_record(gene_id="Rv0006", gene_name="gyrA", annotation=[make_annotation(drug="levofloxacin")])
        ]

        validate_inputs([make_bed_record(locus_tag="Rv0006", gene_name="gyrA")], [], variant_records)

    def test_rejects_a_gene_missing_from_the_gene_database(self, make_variant_record, make_annotation):
        variant_records = [
            make_variant_record(gene_id="RvNotAGene", gene_name="nope", annotation=[make_annotation(drug="levofloxacin")])
        ]

        with pytest.raises(ValueError, match="genes from the results JSON file are missing in the Gene Database: RvNotAGene"):
            validate_inputs([], [], variant_records)

    def test_rejects_an_annotation_drug_missing_from_the_gene_database(self, make_variant_record, make_annotation):
        # gyrA is a real gene, but it is not associated with bedaquiline
        variant_records = [
            make_variant_record(
                gene_id="Rv0006",
                gene_name="gyrA",
                annotation=[make_annotation(drug="bedaquiline")],
                gene_associated_drugs=["levofloxacin"],
            )
        ]

        with pytest.raises(ValueError, match=r"results JSON file are missing in the Gene Database: bedaquiline\|Rv0006\|Rv0006"):
            validate_inputs([], [], variant_records)

    def test_rejects_a_gene_associated_drug_missing_from_the_gene_database(self, make_variant_record, make_annotation):
        variant_records = [
            make_variant_record(
                gene_id="Rv0006",
                gene_name="gyrA",
                annotation=[make_annotation(drug="levofloxacin")],
                gene_associated_drugs=["bedaquiline"],
            )
        ]

        with pytest.raises(ValueError, match=r"results JSON file are missing in the Gene Database: bedaquiline\|Rv0006\|Rv0006"):
            validate_inputs([], [], variant_records)
