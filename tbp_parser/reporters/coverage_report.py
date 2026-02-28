import logging
import pandas as pd
from pathlib import Path

from coverage.coverage_data import TargetCoverage, LocusCoverage
from utils import Configuration, GeneDatabase

logger = logging.getLogger(__name__)


def write_target_coverage_report(
    config: Configuration,
    sample_name: str,
    target_coverage_map: dict[str, TargetCoverage],
    genes_with_valid_deletions: dict[str, list[Variant]],
) -> None:
    """Write a coverage report at the target level.

    Args:
        config: Configuration object
        sample_name: Sample identifier
        target_coverage_map: dict of target_name -> TargetCoverage objects
        genes_with_valid_deletions: dict of gene_name -> list of Variant objects

    Returns:
        None
    """
    rows = []
    for gene_name, coverage in target_coverage_map.items():
        locus_tag = coverage.locus_tag

        qc_warning = (
            f"Contains valid deletion(s): {'; '.join([str(v.nucleotide_change) for v in genes_with_valid_deletions[locus_tag] if coverage.contains_position(v.pos)])}"
            if locus_tag in genes_with_valid_deletions else ""
        )
        rows.append({
            "sample_name": sample_name,
            "locus_tag": coverage.locus_tag,
            "gene_name": gene_name,
            "percent_coverage": f"{coverage.breadth_of_coverage*100:.2f}",
            "average_depth": f"{coverage.average_depth:.2f}",
            "qc_warning": qc_warning,
        })

    df = pd.DataFrame(rows)
    output_path = Path(f"{config.OUTPUT_PREFIX}.target_coverage_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Target coverage report written to {output_path}")

def write_locus_coverage_report(
    config: Configuration,
    sample_name: str,
    locus_coverage_map: dict[str, LocusCoverage],
    genes_with_valid_deletions: dict[str, list[Variant]],
) -> None:
    """Write a coverage report at the locus level.

    Args:
        config: Configuration object
        sample_name: Sample identifier
        locus_coverage_map: dict of locus_tag -> LocusCoverage objects
        genes_with_valid_deletions: dict of gene_name -> list of Variant objects

    Returns:
        None
    """
    rows = []
    for locus_tag, coverage in locus_coverage_map.items():
        qc_warning = (
            f"Contains valid deletion(s): {'; '.join([str(v.nucleotide_change) for v in genes_with_valid_deletions[locus_tag]])}"
            if locus_tag in genes_with_valid_deletions else ""
        )
        rows.append({
            "sample_name": sample_name,
            "locus_tag": locus_tag,
            "gene_name": GeneDatabase.get_gene_name(locus_tag),
            "percent_coverage": f"{coverage.breadth_of_coverage*100:.2f}",
            "average_depth": f"{coverage.average_depth:.2f}",
            "qc_warning": qc_warning,
        })

    df = pd.DataFrame(rows)
    output_path = Path(f"{config.OUTPUT_PREFIX}.locus_coverage_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Locus coverage report written to {output_path}")
