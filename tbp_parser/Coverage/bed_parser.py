import logging
from typing import List
from Coverage import BedRecord

logger = logging.getLogger(__name__)

def parse_bed_file(bed_file: str) -> List[BedRecord]:
    """Parses a BED file and creates BedRecord instances.

    Args:
        bed_file (str): The path to the BED file to parse.
    Returns:
        list[BedRecord]: A list of BedRecord instances parsed from the BED file.
            representing the columns in the BED file.
    """
    bed_records = []
    if not bed_file:
        return bed_records

    logger.debug(f"Parsing BED file: {bed_file}")

    with open(bed_file, 'r') as bf:
        for entry in bf:
            bed_record = BedRecord.from_bed_line(entry)
            bed_records.append(bed_record)

    _validate_unique_bed_records(bed_records)
    logger.debug(f"Parsed {len(bed_records)} records from {bed_file}")
    return bed_records


def _validate_unique_bed_records(bed_records: List[BedRecord]) -> None:
    """Validates that no two BedRecords share the same locus_tag and gene_name.

    Args:
        bed_records: List of BedRecord instances to validate.
    Raises:
        ValueError: If duplicate locus_tag + gene_name combinations are found.
    """
    seen = {}
    duplicates = []

    for record in bed_records:
        key = (record.locus_tag, record.gene_name)
        if key in seen:
            duplicates.append(key)
        else:
            seen[key] = record

    if duplicates:
        raise ValueError(
            f"Duplicate BedRecords found with identical locus_tag and gene_name: "
            f"Records should either be combined into a single entry or split into unique entries."
        )