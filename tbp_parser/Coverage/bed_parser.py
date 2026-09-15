import logging
from typing import List
from tbp_parser.Coverage.bed_record import BedRecord

logger = logging.getLogger(__name__)

def parse_bed_file(bed_file: str, expected_columns: int) -> List[BedRecord]:
    """Parses a BED file and creates BedRecord instances.

    Args:
        bed_file (str): The path to the BED file to parse.
        expected_columns (int): The number of tab-separated columns every row must have.
            `--coverage_bed` and `--err_coverage_bed` need 5 (through gene_name); the
            `build_gene_db --db_bed` file needs 6, because its drug column is the truth set of
            gene/drug associations and `build_gene_database` drops any record without one.
    Returns:
        list[BedRecord]: A list of BedRecord instances parsed from the BED file.
            representing the columns in the BED file.
    Raises:
        ValueError: If any populated row has fewer than `expected_columns` columns.
    """
    bed_records = []
    if not bed_file:
        return bed_records

    logger.debug(f"Parsing BED file: {bed_file}")

    malformed = []
    with open(bed_file, 'r') as bf:
        for number, entry in enumerate(bf, start=1):
            if not entry.strip():
                continue
            if len(entry.strip().split('\t')) < expected_columns:
                malformed.append(number)
                continue
            bed_records.append(BedRecord.from_bed_line(entry))

    if malformed:
        raise ValueError(
            f"{bed_file} requires {expected_columns} tab-separated columns; "
            f"line(s) {', '.join(str(number) for number in malformed)} do not have them"
        )

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