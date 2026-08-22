import argparse
import os
import subprocess
import sys
import logging
import pysam

from tbp_parser.GeneDB.gene_db import GeneDatabase

logger = logging.getLogger(__name__)

def is_file_valid(filename: str) -> str:
    """Checks if an input file is accessible

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    if not os.path.exists(filename) and filename != "-":
        logger.error(f"{filename} cannot be accessed")
        raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    return filename

def is_optional_file_valid(filename: str) -> str:
    """Checks if an optional input file is accessible (no default file provided)

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    if filename != "":
        if not os.path.exists(filename) and filename != "-":
            logger.error(f"{filename} cannot be accessed")
            raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    return filename

def is_bam_index_valid(filename: str) -> str:
    """Checks if there's an associated BAI for the BAM

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    try:
        with pysam.AlignmentFile(filename, "rb") as bam:
            bam.check_index()
    except (OSError, AttributeError) as e:
        logger.error(f"Invalid BAM  for '{filename}': {e}")
        raise argparse.ArgumentTypeError(f"Invalid BAM for '{filename}': {e}")
    except ValueError:
        logger.error("tbp-parser: Generating a BAM index for the input BAM since the BAI appears to be missing / invalid. This could take a while.")
        pysam.index(filename)

    return filename

def is_bed_valid(filename: str) -> str:
    """Checks if the coverage_bed files are accessible

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """

    # check if the necessary columns are present in the BED file -- just count them because we can't really parse it here
    # does this file have at least 5 columns
    if filename != "" and not os.path.exists(filename) and filename != "-":
        logger.error(f"{filename} cannot be accessed")
        raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    else:
        with open(filename, 'r') as bed_file:
            for line in bed_file:
                cols = line.strip().split('\t')
                if len(cols) < 5:
                    logger.error(f"{filename} does not have at least 5 columns as required")
                    raise argparse.ArgumentTypeError("{0} does not have at least 5 columns as required".format(filename))
                break  # only need to check the first line
    return filename

def is_fraction_valid(value: str) -> float:
    """Checks that a percentage-style threshold is a fraction between 0.0 and 1.0

    These thresholds are expressed as fractions (1.0 -> 100%), so a value above 1.0
    can never be met and would silently fail every locus.

    Args:
        value (String): The value to check

    Returns:
        Float: The value as a float if it falls within 0.0 - 1.0
    """
    try:
        fraction = float(value)
    except ValueError:
        logger.error(f"{value} is not a number")
        raise argparse.ArgumentTypeError("{0} is not a number".format(value))

    if not 0.0 <= fraction <= 1.0:
        logger.error(f"{value} must be a fraction between 0.0 and 1.0 (1.0 -> 100%)")
        raise argparse.ArgumentTypeError("{0} must be a fraction between 0.0 and 1.0 (1.0 -> 100%)".format(value))

    return fraction

def is_boundary_valid(boundary_string: str) -> str:
    """Checks if the boundary string for tNGS is valid (two comma-separated numerical values)

    Args:
        boundary_string (String): The boundary string to check
    Returns:
        String: The boundary string if valid
    """
    cols = boundary_string.split(',')
    if len(cols) != 2:
        logger.error(f"{boundary_string} is not formatted correctly; must be two comma-separated values")
        raise argparse.ArgumentTypeError("{0} is not formatted correctly; must be two comma-separated values".format(boundary_string))

    # check if values are numeric
    for val in cols:
        try:
            float(val)
        except ValueError:
            logger.error(f"{boundary_string} is not formatted correctly; both values must be numeric")
            raise argparse.ArgumentTypeError("{0} is not formatted correctly; both values must be numeric".format(boundary_string))

    return boundary_string

def check_bed_for_lims_genes(bed_records, lims_records) -> None:
    """Checks that the provided LIMS format yml file has an associated BedRecord for LIMS report coverage calculations.

    Args:
        bed_records: List of BedRecord objects to check
        lims_records: List of LIMSRecord objects to check
    """

    # It's typical for the lims yaml input file to contain gene names, but the BED file might have irregular
    # gene names representing partial regions. Look for locus tags in the BED file and convert those to gene names to compare with LIMS genes.
    unique_bed_genes = set(record.locus_tag for record in bed_records)
    unique_lims_genes = set()

    missing_from_database = set()
    for rec in lims_records:
        for gene in rec.gene_codes.keys():
            locus_tag = GeneDatabase.get_locus_tag(gene)
            if locus_tag is None:
                missing_from_database.add(gene)
            else:
                unique_lims_genes.add(locus_tag)

    if missing_from_database:
        logger.error(f"The following genes from the LIMS report format yaml file are missing in the Gene Database: {', '.join(missing_from_database)}")
        raise ValueError(f"The following genes from the LIMS report format yaml file are missing in the Gene Database: {', '.join(missing_from_database)}")

    if not unique_lims_genes.issubset(unique_bed_genes):
        missing_locus_tags = unique_lims_genes - unique_bed_genes
        missing_genes_list = [f"{GeneDatabase.get_gene_name(locus_tag)}|{GeneDatabase.get_locus_tag(locus_tag)}" for locus_tag in missing_locus_tags]
        logger.error(f"The following genes from the LIMS report format yaml file are missing in the BED file: {', '.join(missing_genes_list)}")
        raise ValueError(f"The following genes from the LIMS report format yaml file are missing in the BED file: {', '.join(missing_genes_list)}")
    else:
        logger.info("All genes from the LIMS report format yaml file are present in the BED file and Gene Database.")