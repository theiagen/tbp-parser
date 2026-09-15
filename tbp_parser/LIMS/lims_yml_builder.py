import logging
import yaml

from collections import defaultdict

logger = logging.getLogger(__name__)

# TO-DO: Figure out if we can source the drug codes somewhere upstream instead of hardcoding them here.
DRUG_CODES: dict[str, str] = {
  "amikacin": "AMK",
  "bedaquiline": "BDQ",
  "capreomycin": "CAP",
  "clofazimine": "CFZ",
  "cycloserine": "CS",
  "delamanid": "DLM",
  "ethambutol": "EMB",
  "ethionamide": "ETO",
  "isoniazid": "INH",
  "kanamycin": "KAN",
  "levofloxacin": "LFX",
  "linezolid": "LZD",
  "moxifloxacin": "MFX",
  "para-aminosalicylic_acid": "PAS",
  "pretomanid": "PMD",
  "prothionamide": "PTO",
  "pyrazinamide": "PZA",
  "rifampicin": "RIF",
  "rifapentine": "RPT",
  "streptomycin": "STM",
}

def build_lims_fmt(
    gene_database_yml: str,
    output_path: str
) -> None:
    """
    Entry point for the `build_lims_fmt` subcommand.

    Args:
        gene_database_yml (String): The path to the gene database YAML file to derive the format from
        output_path (String): The path of the LIMS report format YAML file to write
    """

    logger.info(f"Building a LIMS report format from `--gene_database_yml` '{gene_database_yml}'")
    lims_report_format = build_lims_report_format(gene_database_yml)
    write_lims_report_format_yml(lims_report_format, output_path)


def _collect_genes_by_drug(gene_database: dict) -> dict[str, set]:
    """
    Inverts a gene database into the set of gene names associated with each drug.

    Args:
        gene_database (dict): A gene database keyed by locus tag

    Returns:
        dict: A map of drugs to their set of associated gene names
    """
    genes_by_drug = defaultdict(set)
    without_drugs = set()

    for locus_tag, entry in gene_database.items():
        gene_name = entry.get("gene_name", locus_tag)
        drugs = entry.get("drugs", [])

        # a gene that is associated with no drug can't appear in a LIMS report
        if not drugs:
            without_drugs.add(gene_name)
            continue

        for drug in drugs:
            genes_by_drug[drug].add(gene_name)

    if without_drugs:
        logger.warning(f"The following {len(without_drugs)} genes have no associated drugs and were left out of the LIMS report format: {', '.join(sorted(without_drugs))}")

    return genes_by_drug


def build_lims_report_format(gene_database_yml: str) -> list[dict]:
    """
    Builds a LIMS report format from a gene database YAML file.

    Args:
        gene_database_yml (String): The path to the gene database YAML file

    Returns:
        list: A LIMS report format, in the same format as a --lims_report_format_yml file
    """
    with open(gene_database_yml, "r") as gene_database_file:
        gene_database = yaml.safe_load(gene_database_file)

    logger.debug(f"Loaded {len(gene_database)} genes from {gene_database_yml}")

    genes_by_drug = _collect_genes_by_drug(gene_database)

    lims_report_format = []
    derived_codes = {}

    for drug in sorted(genes_by_drug):
        drug_code = DRUG_CODES.get(drug)
        if not drug_code:
            drug_code = drug.upper()
            derived_codes[drug] = drug_code

        lims_report_format.append({
            "drug": drug,
            "drug_code": drug_code,
            "gene_codes": {gene: f"{drug_code}_{gene}" for gene in sorted(genes_by_drug[drug], key=str.lower)},
        })

    if derived_codes:
        logger.warning(f"The following {len(derived_codes)} drugs are not in the known drug code map; defulted to their name capitalized: {', '.join(f'{drug} ({code})' for drug, code in sorted(derived_codes.items()))}")

    gene_code_count = sum(len(entry["gene_codes"]) for entry in lims_report_format)
    logger.info(f"Built a LIMS report format with {gene_code_count} gene codes across {len(lims_report_format)} drugs: {', '.join(entry['drug'] for entry in lims_report_format)}")

    return lims_report_format


def write_lims_report_format_yml(
    lims_report_format: list[dict],
    output_path: str
) -> None:
    """
    Writes a LIMS report format to a YAML file.

    Args:
        lims_report_format (list): A LIMS report format, one entry per drug
        output_path (String): The path of the YAML file to write
    """
    with open(output_path, "w") as output_file:
        yaml.dump(lims_report_format, output_file, sort_keys=False, default_flow_style=False)
    logger.info(f"LIMS report format written to '{output_path}'")
