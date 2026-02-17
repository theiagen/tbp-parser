import yaml
import logging
from lims import LIMSRecord, LIMSGeneCode

logger = logging.getLogger(__name__)

def parse_lims_yml_file(lims_report_format_yml: str) -> list[LIMSRecord]:
    """Load the LIMS report format from the YAML file.

    Returns:
        list[LIMSRecord]: List of LIMSRecord objects representing the LIMS report format.
    """
    lims_records = []
    with open(lims_report_format_yml, "r") as f:
        lims_format_list = yaml.safe_load(f)
        for rec in lims_format_list:
            drug=rec["drug"]
            drug_code=rec["drug_code"]

            gene_codes = {}
            for gene, gene_code in rec["gene_codes"].items():
                gene_codes[gene] = LIMSGeneCode(gene_code=gene_code)

            lims_records.append(
                LIMSRecord(
                    drug=drug,
                    drug_code=drug_code,
                    gene_codes=gene_codes
                )
            )

    logger.debug(f"Loaded LIMS report format with {len(lims_records)} drugs")
    return lims_records