import logging
import json

from typing import List, Tuple
from variant import VariantRecord

logger = logging.getLogger(__name__)

def parse_tbprofiler_json(input_json) -> Tuple[List[VariantRecord], str, str, str]:
    """Parses TBProfiler JSON output and returns a list of VariantRecord objects
    which will ultimately be converted to Variant objects.

    Returns:
        Tuple[List[VariantRecord], str, str, str]: A tuple containing
            - list of VariantRecord objects parsed from the TBProfiler JSON output
            - the sample ID
            - the lineage ID
            - the sublineage ID
    """
    with open(input_json) as f:
        json_data = json.load(f)

    SAMPLE_ID = json_data.get("id")
    LINEAGE_ID = json_data.get("main_lineage")
    SUBLINEAGE_ID = json_data.get("sub_lineage")

    logger.debug(f"Parsing TBProfiler JSON: {input_json}")
    logger.debug(f"Parsing `dr_variants` and `other_variants` JSON entries for sample ID: {SAMPLE_ID}")
    logger.debug(f"Found {len(json_data.get('dr_variants', []))} `dr_variants` and {len(json_data.get('other_variants', []))} `other_variants`")

    # process both dr_variants and other_variants
    all_variant_records = []
    for variant_record in json_data.get("dr_variants", []) + json_data.get("other_variants", []):
        variant_record["sample_id"] = SAMPLE_ID
        variant_record = VariantRecord(**variant_record)
        all_variant_records.append(
            variant_record
        )
    return all_variant_records, SAMPLE_ID, LINEAGE_ID, SUBLINEAGE_ID