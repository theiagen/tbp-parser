import logging
import pandas as pd
from datetime import datetime
from pathlib import Path

from Variant import Variant
from Utilities import Configuration, apply_find_and_replace

logger = logging.getLogger(__name__)

# Looker preserves interim designations
LOOKER_RESISTANCE_RANKING = {
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
    variants: list[Variant],
    lims_lineage: str,
    sample_id: str,
    detected_lineage: str,
) -> None:
    """Write the Looker report from processed Variant objects.

    Args:
        variants: All processed Variant objects
        lims_lineage: Lineage string from LIMS
        sample_id: Sample ID string
        detected_lineage: Detected lineage string from TBProfiler

    Returns:
        None
    """
    config = Configuration.get_instance()

    drug_map = {}
    # If there are multiple variants for the same drug, we want to report the highest ranked interpretation
    for variant in variants:
        if variant.drug not in drug_map:
            drug_map[variant.drug] = variant.looker_interpretation
        else:
            current_rank = LOOKER_RESISTANCE_RANKING.get(drug_map[variant.drug], -1)
            new_rank = LOOKER_RESISTANCE_RANKING.get(variant.looker_interpretation or "NA", -1) # "NA" is redundant but pylance was complaining
            if new_rank > current_rank:
                drug_map[variant.drug] = variant.looker_interpretation

    data = {
        "sample_id": sample_id,
        "output_seq_method_type": config.SEQUENCING_METHOD,
    }

    data = {
        **data,
        **{
            k: drug_map[k]
            for k in sorted(drug_map)
            if k not in ["cycloserine", "delamanid", "para-aminosalicylic_acid"] # remove these drugs from the report for now, to be consistent with old report
        }
    }

    data["lineage"] = detected_lineage
    data["ID"] = lims_lineage
    data["analysis_date"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    data["operator"] = config.OPERATOR

    df = pd.DataFrame([data])
    df = apply_find_and_replace(df, config.FIND_AND_REPLACE)
    output_path = Path(f"{config.OUTPUT_PREFIX}.looker_report.csv")
    df.to_csv(output_path, index=False)

    logger.info(f"Looker report written to {output_path}")
