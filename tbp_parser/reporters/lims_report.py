"""LIMS report orchestrator — handles I/O, iteration, and CSV writing.

All decision logic lives in lims_rules.py. This module loads configuration,
filters variants, and calls lims_rules functions to populate the report.
"""

import logging
import pandas as pd
from datetime import datetime
from pathlib import Path

from lims import LIMSRecord
from utils import Configuration

logger = logging.getLogger(__name__)

def write_lims_report(
    config: Configuration,
    lims_records: list[LIMSRecord],
    lims_lineage: str,
    sample_id: str,
    detected_lineage: str,
) -> None:
    """Write the LIMS report from processed Variant objects.

    Args:
        config: Configuration object
        lims_records: List of LIMSRecord objects
        lims_lineage: Lineage information for the LIMS report
        sample_id: Sample identifier
        detected_lineage: Detected lineage identifier

    Returns:
        None
    """
    # Build the report as a dictionary first, then convert to DataFrame for writing
    lims_report = {
        "MDL sample accession numbers": sample_id,
        "M_DST_A01_ID": lims_lineage,
    }

    for record in lims_records:
        lims_report[record.drug_code] = record.drug_target_value if record.drug_target_value is not None else ""
        for gene, gene_code in record.gene_codes.items():
            lims_report[gene_code.gene_code] = gene_code.gene_target_value if gene_code.gene_target_value is not None else ""

    # Add metadata
    lims_report["Analysis date"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    lims_report["Operator"] = config.OPERATOR
    lims_report["M_DST_O01_Lineage"] = detected_lineage

    # Write CSV and transposed CSV
    df = pd.DataFrame([lims_report])
    output_path = Path(f"{config.OUTPUT_PREFIX}.lims_report.csv")
    df.to_csv(output_path, index=False)

    transposed_path = Path(f"{config.OUTPUT_PREFIX}.lims_report.transposed.csv")
    df.T.to_csv(transposed_path, header=False)

    logger.info(f"LIMS report written to {output_path}")
