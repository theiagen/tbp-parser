import logging
import pandas as pd
from datetime import datetime
from pathlib import Path

from variant import Variant
from utils.config import Configuration

logger = logging.getLogger(__name__)

# Looker preserves interim designations
RESISTANCE_RANKING = {
    "R": 6,
    "R-Interim": 5,
    "U": 4,
    "S-Interim": 3,
    "S": 2,
    "WT": 1,
    "Insufficient Coverage": 0,
    "NA": -1,
}


def write_looker_report(
    config: Configuration,
    variants: list[Variant],
    low_depth_genes: list[str],
    genes_with_valid_deletions: list[str],
    lineage: str,
    lineage_english: str,
) -> Path:
    """Write the Looker report from processed Variant objects.

    Args:
        config: Configuration object
        variants: All processed Variant objects
        low_depth_genes: Genes that failed locus coverage QC
        genes_with_valid_deletions: Genes with QC-passing deletions
        lineage: Raw lineage string from TBProfiler
        lineage_english: Human-readable lineage string

    Returns:
        Path to the written CSV file
    """
    lims_format = config.load_lims_report_format()

    # Build reverse mapping: drug -> list of gene names
    drugs_to_genes = {}
    for drug, drug_info in lims_format.items():
        for _code, genes in drug_info.items():
            for gene_name in genes:
                if drug not in drugs_to_genes:
                    drugs_to_genes[drug] = []
                drugs_to_genes[drug].append(gene_name)

    antimicrobial_resistances = {}
    for antimicrobial, genes in drugs_to_genes.items():
        drug_variants = [v for v in variants if v.drug == antimicrobial]
        looker_interpretations = [
            v.looker_interpretation for v in drug_variants
            if v.looker_interpretation
        ]

        try:
            max_looker = max(
                looker_interpretations,
                key=lambda x: RESISTANCE_RANKING.get(x, -1)
            )
        except ValueError:
            max_looker = "NA"

        antimicrobial_resistances[antimicrobial] = max_looker

        # Override with insufficient coverage if any gene has low depth (and non-R)
        if max_looker != "R":
            for gene in genes:
                if gene in low_depth_genes and gene not in genes_with_valid_deletions:
                    antimicrobial_resistances[antimicrobial] = "Insufficient coverage in locus"
                    break

    sample_name = variants[0].sample_id if variants else ""
    data = {
        "sample_id": sample_name,
        "output_seq_method_type": config.SEQUENCING_METHOD,
    }
    for antimicrobial in sorted(antimicrobial_resistances.keys()):
        data[antimicrobial] = antimicrobial_resistances[antimicrobial]

    data["lineage"] = lineage
    data["ID"] = lineage_english
    data["analysis_date"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    data["operator"] = config.OPERATOR

    df = pd.DataFrame([data])
    output_path = Path(f"{config.OUTPUT_PREFIX}.looker_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Looker report written to {output_path}")
    return output_path
