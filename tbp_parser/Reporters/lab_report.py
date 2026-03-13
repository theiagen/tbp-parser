import logging
from typing import Optional
import pandas as pd
from pathlib import Path

from tbp_parser.Variant.variant import Variant
from tbp_parser.Utilities.config import Configuration, apply_find_and_replace

logger = logging.getLogger(__name__)

# Define the columns and their corresponding Variant attributes
COLUMNS = {
    "sample_id": "sample_id",
    "tbprofiler_gene_name": "gene_name",
    "tbprofiler_locus_tag": "gene_id",
    "tbprofiler_variant_substitution_type": "type",
    "tbprofiler_variant_substitution_nt": "nucleotide_change",
    "tbprofiler_variant_substitution_aa": "protein_change",
    "genomic_start_pos": "absolute_start",
    "genomic_end_pos": "absolute_end",
    "confidence": "confidence",
    "antimicrobial": "drug",
    "looker_interpretation": "looker_interpretation",
    "mdl_interpretation": "mdl_interpretation",
    "depth": "depth",
    "frequency": "freq",
    "read_support": "read_support",
    "rationale": "rationale",
    "warning": "warning",
    "gene_tier": "gene_tier",
    "source": "source",
    "tbdb_comment": "comment",
}

def get_report_fmt(variant: Variant, attribute: str) -> Optional[str]:
    """Format attribute for report output, converting missing/special attributes to specific values.
    Can be used for standardizing unreported variants (WT/NA) in reports.

    Args:
        attribute: The attribute name to format.

    Returns:
        The formatted attribute value as a string.
    """
    value = getattr(variant, attribute)
    if attribute == "warning":
        return "; ".join(sorted(value))
    else:
        return str(value)

def sort_lab_df(df: pd.DataFrame) -> pd.DataFrame:
    """Sort the laboratorian report DataFrame.

    Args:
        df: DataFrame of the laboratorian report

    Returns:
        Sorted DataFrame
    """
    # split into three DataFrames based on nucleotide substitution -> other, WT, NA
    df_other = df[~df["tbprofiler_variant_substitution_nt"].isin(["WT", "NA"])]
    df_wt = df[df["tbprofiler_variant_substitution_nt"] == "WT"]
    df_na = df[df["tbprofiler_variant_substitution_nt"] == "NA"]

    # then sort each DataFrame by gene name (case insensitive)
    df_other = df_other.sort_values(by=["tbprofiler_gene_name", "tbprofiler_variant_substitution_nt", "antimicrobial"], key=lambda col: col.str.lower()) # type: ignore
    df_wt = df_wt.sort_values(by=["tbprofiler_gene_name", "tbprofiler_variant_substitution_nt", "antimicrobial"], key=lambda col: col.str.lower()) # type: ignore
    df_na = df_na.sort_values(by=["tbprofiler_gene_name", "tbprofiler_variant_substitution_nt", "antimicrobial"], key=lambda col: col.str.lower()) # type: ignore

    # concatenate them back together
    df_all = pd.concat([df_other, df_wt, df_na], ignore_index=True)
    return df_all

def write_laboratorian_report(variants: list[Variant]) -> None:
    """Write the laboratorian report from processed Variant objects.

    Args:
        variants: List of Variant objects with interpretation and QC applied

    Returns:
        None
    """
    config = Configuration.get_instance()
    rows = []
    for variant in variants:
        row = {}
        for col, attr in COLUMNS.items():
            row[col] = get_report_fmt(variant, attr)

        rows.append(row)

    df = pd.DataFrame(rows, columns=list(COLUMNS.keys()))
    df = apply_find_and_replace(df, config.FIND_AND_REPLACE)
    df = sort_lab_df(df)
    output_path = Path(f"{config.OUTPUT_PREFIX}.laboratorian_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Laboratorian report written to {output_path}")