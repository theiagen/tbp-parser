import logging
from typing import Optional
import pandas as pd
from pathlib import Path

from variant import Variant
from utils.config import Configuration

logger = logging.getLogger(__name__)

# Define the columns and their corresponding Variant attributes
COLUMNS = {
    "sample_id": "sample_id",
    "tbprofiler_gene_name": "gene_name",
    "tbprofiler_locus_tag": "gene_id",
    "tbprofiler_variant_substitution_type": "type",
    "tbprofiler_variant_substitution_nt": "nucleotide_change",
    "tbprofiler_variant_substitution_aa": "protein_change",
    # "tbprofiler_variant_position": "pos",
    # "variant_genomic_start_pos": "absolute_start",
    # "variant_genomic_end_pos": "absolute_end",
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
    if not value:
        # conditions for reporting "NA"
        if attribute == "protein_change":
            return "NA"
        # conditions for reporting empty strings ""
        elif attribute in ["source", "comment", "rationale"]:
            return ""
    # special reporting for non-empty attributes
    elif attribute == "warning":
        return "; ".join(sorted(value))
    # report negative numeric values as "NA" (used for depth, pos, freq, read_support in `unreported_variants``)
    elif isinstance(value, (int, float)) and value < 0:
        return "NA"
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
    df_other = df_other.sort_values(by="tbprofiler_gene_name", key=lambda col: col.str.lower()) # type: ignore
    df_wt = df_wt.sort_values(by="tbprofiler_gene_name", key=lambda col: col.str.lower()) # type: ignore
    df_na = df_na.sort_values(by="tbprofiler_gene_name", key=lambda col: col.str.lower()) # type: ignore

    # concatenate them back together
    df_all = pd.concat([df_other, df_wt, df_na], ignore_index=True)
    return df_all

def write_laboratorian_report(config: Configuration, variants: list[Variant]) -> Path:
    """Write the laboratorian report from processed Variant objects.

    Args:
        config: Configuration object (used for OUTPUT_PREFIX)
        variants: List of Variant objects with interpretation and QC applied

    Returns:
        Path to the written CSV file
    """
    rows = []
    for variant in variants:
        row = {}
        for col, attr in COLUMNS.items():
            row[col] = get_report_fmt(variant, attr)

            # remove this!!! just for comparing to old reports right now
            if row[col] == "mmpR5":
                row[col] = "Rv0678"

        rows.append(row)

    df = pd.DataFrame(rows, columns=COLUMNS.keys())
    df = sort_lab_df(df)
    output_path = Path(f"{config.OUTPUT_PREFIX}.laboratorian_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Laboratorian report written to {output_path}")
    return output_path