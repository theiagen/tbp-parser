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
    for gene_name, target_coverage in target_coverage_map.items():

        # if ERR coverage exists, only report deletions if they are within the ERR region
        deletions = target_coverage.err_coverage.valid_deletions if target_coverage.err_coverage else target_coverage.valid_deletions
        qc_warning = f"Contains valid deletion(s): {'; '.join(str(v.nucleotide_change) for v in deletions)}" if deletions else ""

        # Even if ERR coverage exists, still report the breadth and depth for the full target region (not just the ERR)
        rows.append({
            "sample_name": sample_name,
            "locus_tag": target_coverage.locus_tag,
            "gene_name": gene_name,
            "percent_coverage": f"{target_coverage.breadth_of_coverage*100:.3f}",
            "average_depth": f"{target_coverage.average_depth:.3f}",
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
    for locus_tag, locus_coverage in locus_coverage_map.items():

        # if ERR coverage exists, only report deletions if they are within the ERR region
        deletions = locus_coverage.err_coverage.valid_deletions if locus_coverage.err_coverage else locus_coverage.valid_deletions
        qc_warning = f"Contains valid deletion(s): {'; '.join(str(v.nucleotide_change) for v in deletions)}" if deletions else ""

        rows.append({
            "sample_name": sample_name,
            "locus_tag": locus_tag,
            "gene_name": GeneDatabase.get_gene_name(locus_tag),
            "percent_coverage": f"{locus_coverage.breadth_of_coverage*100:.3f}",
            "average_depth": f"{locus_coverage.average_depth:.3f}",
            "qc_warning": qc_warning,
        })

    df = pd.DataFrame(rows)
    output_path = Path(f"{config.OUTPUT_PREFIX}.locus_coverage_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Locus coverage report written to {output_path}")
