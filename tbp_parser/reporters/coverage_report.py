import logging
import pandas as pd
from pathlib import Path

from coverage import GeneCoverage, LocusCoverage
from utils.config import Configuration

logger = logging.getLogger(__name__)


def write_gene_coverage_report(
    config: Configuration,
    sample_name: str,
    gene_coverage_map: dict[str, GeneCoverage],
) -> Path:
    """Write a coverage report at the gene level.

    Args:
        config: Configuration object
        sample_name: Sample identifier
        gene_coverage_map: dict of gene_name -> GeneCoverage objects

    Returns:
        Path to the written CSV file
    """
    rows = []
    for gene_name, coverage in gene_coverage_map.items():
        rows.append({
            "sample_name": sample_name,
            "gene_name": gene_name,
            "percent_coverage": coverage.breadth_of_coverage,
            "average_depth": coverage.average_depth,
        })

    df = pd.DataFrame(rows)
    output_path = Path(f"{config.OUTPUT_PREFIX}.gene_coverage_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Gene coverage report written to {output_path}")
    return output_path


def write_locus_coverage_report(
    config: Configuration,
    sample_name: str,
    locus_coverage_map: dict[str, LocusCoverage],
) -> Path:
    """Write a coverage report at the locus level.

    Args:
        config: Configuration object
        sample_name: Sample identifier
        locus_coverage_map: dict of locus_tag -> LocusCoverage objects

    Returns:
        Path to the written CSV file
    """
    rows = []
    for locus_tag, coverage in locus_coverage_map.items():
        rows.append({
            "sample_name": sample_name,
            "locus_tag": locus_tag,
            "percent_coverage": coverage.breadth_of_coverage,
            "average_depth": coverage.average_depth,
        })

    df = pd.DataFrame(rows)
    output_path = Path(f"{config.OUTPUT_PREFIX}.locus_coverage_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Locus coverage report written to {output_path}")
    return output_path
