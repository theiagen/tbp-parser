import yaml

from tbp_parser.GeneDB.gene_db_builder import (
    build_gene_database,
    write_gene_database_yml,
)


class TestBuildGeneDatabase:

    def test_builds_an_entry_from_a_bed_record(self, make_bed_record):
        bed_records = [make_bed_record(locus_tag="Rv0006", gene_name="gyrA", drugs=["levofloxacin", "moxifloxacin"])]

        gene_database = build_gene_database(bed_records)

        assert gene_database == {
            "Rv0006": {
                "locus_tag": "Rv0006",
                "gene_name": "gyrA",
                "tier": "Tier 1",
                "promoter_region": [-35, -1],
                "drugs": ["levofloxacin", "moxifloxacin"],
            }
        }

    def test_skips_records_without_drugs(self, make_bed_record):
        # the drug column is the truth set of gene/drug interactions, so a region without one
        # cannot be reported on and does not become an entry.
        bed_records = [
            make_bed_record(locus_tag="Rv0006", gene_name="gyrA", drugs=["levofloxacin"]),
            make_bed_record(locus_tag="Rv0001", gene_name="dnaA"),
        ]

        gene_database = build_gene_database(bed_records)

        assert list(gene_database) == ["Rv0006"]

    def test_carries_aliases_for_the_rrna_genes(self, make_bed_record):
        bed_records = [make_bed_record(locus_tag="EBG00000313325", gene_name="rrs", drugs=["amikacin"])]

        gene_database = build_gene_database(bed_records)

        assert gene_database["EBG00000313325"]["aliases"] == ["Rvnr01", "MTB000019"]

    def test_defaults_the_metadata_of_an_unknown_gene(self, make_bed_record):
        # a gene the metadata does not name still gets an entry, so it stays reportable
        bed_records = [make_bed_record(locus_tag="Rv9999", gene_name="mysteryGene", drugs=["isoniazid"])]

        gene_database = build_gene_database(bed_records)

        assert gene_database["Rv9999"] == {
            "locus_tag": "Rv9999",
            "gene_name": "mysteryGene",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["isoniazid"],
        }

    def test_sorts_and_deduplicates_drugs_within_an_entry(self, make_bed_record):
        bed_records = [make_bed_record(locus_tag="Rv0667", gene_name="rpoB", drugs=["rifapentine", "rifampicin", "rifampicin"])]

        gene_database = build_gene_database(bed_records)

        assert gene_database["Rv0667"]["drugs"] == ["rifampicin", "rifapentine"]


class TestWriteGeneDatabaseYml:

    def test_round_trips_through_the_gene_database(self, tmp_path, make_bed_record):
        bed_records = [
            make_bed_record(locus_tag="Rv0006", gene_name="gyrA", drugs=["levofloxacin"]),
            make_bed_record(locus_tag="EBG00000313325", gene_name="rrs", drugs=["amikacin"]),
        ]
        output_path = tmp_path / "gene-database.yml"

        gene_database = build_gene_database(bed_records)
        write_gene_database_yml(gene_database, str(output_path))

        assert yaml.safe_load(output_path.read_text()) == gene_database

    def test_sorts_entries_by_locus_tag(self, tmp_path, make_bed_record):
        # entries are sorted so that diffs between database versions stay readable.
        bed_records = [
            make_bed_record(locus_tag="Rv0006", gene_name="gyrA", drugs=["levofloxacin"]),
            make_bed_record(locus_tag="EBG00000313325", gene_name="rrs", drugs=["amikacin"]),
        ]
        output_path = tmp_path / "gene-database.yml"

        gene_database = build_gene_database(bed_records)
        write_gene_database_yml(gene_database, str(output_path))

        # the records are deliberately passed in unsorted order, so this fails if the sort is dropped
        assert list(gene_database) == ["Rv0006", "EBG00000313325"]
        assert list(yaml.safe_load(output_path.read_text())) == ["EBG00000313325", "Rv0006"]

    def test_preserves_key_order_within_an_entry(self, tmp_path, make_bed_record):
        # sort_keys=False keeps each entry readable in the written file, locus_tag first
        bed_records = [make_bed_record(locus_tag="Rv0006", gene_name="gyrA", drugs=["levofloxacin"])]
        output_path = tmp_path / "gene-database.yml"

        gene_database = build_gene_database(bed_records)
        write_gene_database_yml(gene_database, str(output_path))

        written = yaml.safe_load(output_path.read_text())
        assert list(written["Rv0006"]) == ["locus_tag", "gene_name", "tier", "promoter_region", "drugs"]