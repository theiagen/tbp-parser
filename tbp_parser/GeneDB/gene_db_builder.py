from __future__ import annotations
import logging
import yaml

from typing import TYPE_CHECKING
from tbp_parser.GeneDB.gene_db_metadata import GENE_DB_METADATA

if TYPE_CHECKING:
    from tbp_parser.Coverage.bed_record import BedRecord

logger = logging.getLogger(__name__)

def build_gene_db(
    bed_records: list[BedRecord],
    output_path: str
) -> None:
    """
    Entry point for the `build_gene_db` subcommand.

    Args:
        bed_records (list[BedRecord]): The records parsed from the `--db_bed` file
        output_path (String): The path of the gene database YAML file to write
    """

    logger.info(f"Building a gene database from {len(bed_records)} `--db_bed` BedRecords")
    gene_database = build_gene_database(bed_records)
    write_gene_database_yml(gene_database, f"{output_path}")


def build_gene_database(bed_records: list[BedRecord]) -> dict:
    """
    Builds a gene database from the `--db_bed` file the variants were called against.

    The `--db_bed` file is the truth set for which genes exist and which drugs each is associated
    with, so its drug column is taken verbatim.

    Args:
        bed_records (list[BedRecord]): The records parsed from the `--db_bed` file

    Returns:
        dict: A gene database keyed by locus tag, in the same format as a --gene_database_yml file
    """
    locus_tag2gene_name = {record.locus_tag: record.gene_name for record in bed_records}
    locus_tag2drugs = {record.locus_tag: record.drugs for record in bed_records if record.drugs}

    gene_database = {}
    supplemented = set()

    for locus_tag, drugs in locus_tag2drugs.items():
        metadata = GENE_DB_METADATA.get(locus_tag, {})
        if metadata:
            supplemented.add(locus_tag)

        entry = {
            "locus_tag": locus_tag,
            "gene_name": locus_tag2gene_name.get(locus_tag, locus_tag),
            "tier": metadata.get("tier", "NA"),
            "promoter_region": metadata.get("promoter_region", []),
            "drugs": sorted(set(drugs)),
        }
        if metadata.get("aliases"):
            entry["aliases"] = list(metadata["aliases"])

        gene_database[locus_tag] = entry

    if supplemented:
        logger.warning(f"A total of {len(supplemented)} genes were supplemented with `tier` and `promoter region` metadata")

    all_drugs = sorted({drug for entry in gene_database.values() for drug in entry["drugs"]})
    logger.info(f"Built a gene database with {len(gene_database)} genes across {len(all_drugs)} drugs: {', '.join(all_drugs)}")

    return gene_database


def write_gene_database_yml(
    gene_database: dict,
    output_path: str
) -> None:
    """
    Writes a gene database to a YAML file.

    Entries are sorted by locus tag so that diffs between database versions stay readable, while
    `sort_keys=False` preserves the key order within each entry.

    Args:
        gene_database (dict): A gene database keyed by locus tag
        output_path (String): The path of the YAML file to write
    """
    with open(output_path, "w") as output_file:
        yaml.dump(dict(sorted(gene_database.items())), output_file, sort_keys=False, default_flow_style=None)
    logger.info(f"Gene database written to '{output_path}'")
    return