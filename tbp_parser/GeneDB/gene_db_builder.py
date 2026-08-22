import argparse
import json
import logging
import yaml

from collections import defaultdict
from tbp_parser.GeneDB.gene_db_metadata import GENE_DB_METADATA

logger = logging.getLogger(__name__)

def build_gene_db(options: argparse.Namespace) -> None:
    """
    Entry point for the `build_gene_db` subcommand.

    Args:
        options (argparse.Namespace): an object with the input arguments provided at runtime
    """

    logger.info(f"Building a gene database from --db '{options.db}' and --input_json '{options.input_json}'")
    gene_database = build_gene_database(options.db, options.input_json)
    write_gene_database_yml(gene_database, f"{options.output}")


def _load_gene_name2locus_tag(input_json_path: str) -> dict[str, str]:
    """
    Loads the `gene_name2locus_tag` field from a TBProfiler (v6.7.0) results JSON file.

    Args:
        input_json_path (String): The path to the TBProfiler results JSON file

    Returns:
        dict: A map of gene name to locus tag
    """
    with open(input_json_path) as input_json:
        results = json.load(input_json)

    gene_name2locus_tag = results.get("gene_name2locus_tag")
    if not gene_name2locus_tag:
        logger.error(f"{input_json_path} does not contain a `gene_name2locus_tag` map; this map is only emitted by TBProfiler v6.7.0 and later")
        raise ValueError(f"{input_json_path} does not contain a `gene_name2locus_tag` map; this map is only emitted by TBProfiler v6.7.0 and later")

    logger.debug(f"Loaded {len(gene_name2locus_tag)} gene name to locus tag mappings from {input_json_path}")
    return gene_name2locus_tag


def _collect_drugs_by_locus_tag(
    db_json_path: str,
    gene_name2locus_tag: dict[str, str]
) -> dict[str, set]:
    """
    Collects the set of drugs associated with each locus tag in the TBProfiler database.

    The `mutations.json` database is the truth set for which genes exist, so every one of its entries is kept.
    `gene_name2locus_tag` is used only to normalize the gene name field if unspecified in the `mutations.json` database.

    Args:
        db_json_path (String): The path to the TBProfiler database JSON file
        gene_name2locus_tag (dict): A map of gene name to locus tag

    Returns:
        dict: A map of locus tag to its set of associated drugs
    """
    with open(db_json_path) as db_json:
        database = json.load(db_json)

    locus_tag2drugs = defaultdict(set)

    for gene, mutations in database.items():
        # the database is normally keyed by locus tag but this also checks for gene name
        locus_tag = gene_name2locus_tag.get(gene, gene)

        for variant, annotations in mutations.items():
            for annotation in annotations.get("annotations", []):
                drug = annotation.get("drug")
                # only creates an entry if drug exists for that locus_tag
                if drug:
                    locus_tag2drugs[locus_tag].add(drug)
                else:
                    logger.warning(f"No drugs associated with {locus_tag} {variant}")

    return locus_tag2drugs


def build_gene_database(
    db_json_path: str,
    input_json_path: str
) -> dict:
    """
    Builds a gene database from the TBProfiler database and results JSON files used to call variants.

    Args:
        db_json_path (String): The path to the TBProfiler database JSON file
        input_json_path (String): The path to the TBProfiler results JSON file

    Returns:
        dict: A gene database keyed by locus tag, in the same format as a --gene_database_yml file
    """
    gene_name2locus_tag = _load_gene_name2locus_tag(input_json_path)
    locus_tag2gene_name = {locus_tag: gene_name for gene_name, locus_tag in gene_name2locus_tag.items()}

    locus_tag2drugs = _collect_drugs_by_locus_tag(db_json_path, gene_name2locus_tag)

    gene_database = {}
    defaulted = set()

    for locus_tag, drugs in locus_tag2drugs.items():
        metadata = GENE_DB_METADATA.get(locus_tag, {})
        if metadata:
            defaulted.add(locus_tag)

        entry = {
            "locus_tag": locus_tag,
            "gene_name": locus_tag2gene_name.get(locus_tag, locus_tag),
            "tier": metadata.get("tier", "NA"),
            "promoter_region": metadata.get("promoter_region", []),
            "drugs": sorted(list(drugs)),
        }
        if metadata.get("aliases"):
            entry["aliases"] = list(metadata["aliases"])

        gene_database[locus_tag] = entry

    if defaulted:
        logger.warning(f"A total of {len(defaulted)} genes are given default `tier` and `promoter region` metadata")

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
