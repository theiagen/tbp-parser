import logging
import pandas as pd
from pathlib import Path

from coverage import TargetCoverage, LocusCoverage
from utils import Configuration

logger = logging.getLogger(__name__)


def write_target_coverage_report(
    config: Configuration,
    sample_name: str,
    target_coverage_map: dict[str, TargetCoverage],
) -> None:
    """Write a coverage report at the target level.

    Args:
        config: Configuration object
        sample_name: Sample identifier
        target_coverage_map: dict of target_name -> TargetCoverage objects

    Returns:
        None
    """
    rows = []
    for gene_name, coverage in target_coverage_map.items():
        rows.append({
            "sample_name": sample_name,
            "gene_name": gene_name,
            "percent_coverage": coverage.breadth_of_coverage,
            "average_depth": coverage.average_depth,
        })

    df = pd.DataFrame(rows)
    output_path = Path(f"{config.OUTPUT_PREFIX}.target_coverage_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Target coverage report written to {output_path}")

def write_locus_coverage_report(
    config: Configuration,
    sample_name: str,
    locus_coverage_map: dict[str, LocusCoverage],
) -> None:
    """Write a coverage report at the locus level.

    Args:
        config: Configuration object
        sample_name: Sample identifier
        locus_coverage_map: dict of locus_tag -> LocusCoverage objects

    Returns:
        None
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
