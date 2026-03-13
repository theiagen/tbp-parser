import logging
import pandas as pd
from pathlib import Path

from Coverage import TargetCoverage, LocusCoverage
from Utilities import Configuration, GeneDatabase, apply_find_and_replace

logger = logging.getLogger(__name__)


def write_coverage_report(
    sample_name: str,
    coverage_map: dict[str, TargetCoverage] | dict[str, LocusCoverage],
) -> None:
    """Write a coverage report at the target level.

    Args:
        sample_name: Sample identifier
        coverage_map: either dict of gene_name/locus_tag -> TargetCoverage/LocusCoverage objects

    Returns:
        None
    """
    config = Configuration.get_instance()

    rows = []
    for key, coverage in coverage_map.items():
        if isinstance(coverage, LocusCoverage):
            locus_tag = key
            gene_name = GeneDatabase.get_gene_name(key)
            output_path = Path(f"{config.OUTPUT_PREFIX}.locus_coverage_report.csv")

        else: # TargetCoverage
            locus_tag = coverage.locus_tag
            gene_name = key
            output_path = Path(f"{config.OUTPUT_PREFIX}.target_coverage_report.csv")

        # if ERR coverage exists, only report deletions if they are within the ERR region
        deletions = coverage.err_coverage.valid_deletions if coverage.err_coverage else coverage.valid_deletions
        if deletions:
            qc_warning = (
                f"Contains valid deletion(s) in " +
                f"{'ERR ' if coverage.err_coverage else 'full '}" +
                f"{'locus region: ' if isinstance(coverage, LocusCoverage) else 'target region: '}" +
                f"{'; '.join(str(v.nucleotide_change) for v in deletions)}"
            )
        else:
            qc_warning = ""

        # Even if ERR coverage exists, still report the breadth and depth for the full target region (not just the ERR)
        rows.append({
            "sample_name": sample_name,
            "locus_tag": locus_tag,
            "gene_name": gene_name,
            "percent_coverage": f"{coverage.breadth_of_coverage*100:.3f}",
            "average_depth": f"{coverage.average_depth:.3f}",
            **({"err_percent_coverage": f"{coverage.err_coverage.breadth_of_coverage*100:.3f}"} if coverage.err_coverage else {}),
            **({"err_average_depth": f"{coverage.err_coverage.average_depth:.3f}"} if coverage.err_coverage else {}),
            "qc_warning": qc_warning,
        })

    df = pd.DataFrame(rows)
    df = apply_find_and_replace(df, config.FIND_AND_REPLACE)
    df.to_csv(output_path, index=False)  # type: ignore

    logger.info(f"Coverage report written to {output_path}") # type: ignore
