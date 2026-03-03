import argparse
import os
import subprocess
import sys
import logging

from Utilities import GeneDatabase

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

def is_bam_valid(filename: str) -> str:
    """Checks if the input BAM is valid and accessible and if there is an associated BAI

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    if not os.path.exists(filename) and (filename != "-" and not filename.endswith(".bam")):
        logger.error(f"{filename} cannot be accessed or is missing the BAM extension")
        raise argparse.ArgumentTypeError("{0} cannot be accessed or is missing the BAM extension".format(filename))
    bai_filename = filename + ".bai"
    if not os.path.exists(bai_filename) and (bai_filename != "-"):
        logger.error(f"Cannot find the BAI file for this BAM: {bai_filename}")
        raise argparse.ArgumentTypeError("Cannot find the BAI file for this BAM")
    return filename

def is_bed_valid(filename: str) -> str:
    """Checks if the TBDB BED file is accessible

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """

    # check if the necessary columns are present in the BED file -- just count them because we can't really parse it here
    # does this file have at least 6 columns
    if filename != "" and not os.path.exists(filename) and filename != "-":
        logger.error(f"{filename} cannot be accessed")
        raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    else:
        with open(filename, 'r') as bed_file:
            for line in bed_file:
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    logger.error(f"{filename} does not have at least 6 columns as required")
                    raise argparse.ArgumentTypeError("{0} does not have at least 6 columns as required".format(filename))
                break  # only need to check the first line
    return filename

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

def check_dependency_exists() -> None:
    """This function confirms that samtools is installed and available; if it is
    not, the program exits with an error message.
    """
    print("Checking for samtools installation...")
    try:
      result = subprocess.run(
          ["samtools", "--version"],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE,
          encoding='utf-8',
          errors='replace',
          check=True,
          text=True,
      )
      version = result.stdout.split('\n')[0]
      logger.info(f"Found samtools: {version}")
    except FileNotFoundError:
        # samtools not found
        logger.error("Error: samtools not found. Please install samtools and try again.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        # samtools found but failed to execute properly
        logger.error(f"samtools found but failed to run properly: {e.stderr.strip()}")
        sys.exit(1)

def check_bed_for_lims_genes(bed_records, lims_records) -> None:
    """Checks that the provided LIMS format yml file has an associated BedRecord for LIMS report coverage calculations.

    Args:
        bed_records: List of BedRecord objects to check
        lims_records: List of LIMSRecord objects to check
    """

    # It's typical for the lims yaml input file to contain gene names, but the BED file might have irregular
    # gene names representing partial regions. Look for locus tags in the BED file and convert those to gene names to compare with LIMS genes.
    unique_bed_genes = set(record.locus_tag for record in bed_records)
    unique_lims_genes = set([
        GeneDatabase.get_locus_tag(gene)
        for rec in lims_records
        for gene in rec.gene_codes.keys()
    ])

    if not unique_lims_genes.issubset(unique_bed_genes):
        missing_locus_tags = unique_lims_genes - unique_bed_genes
        missing_genes_list = [f"{GeneDatabase.get_gene_name(locus_tag)}|{GeneDatabase.get_locus_tag(locus_tag)}" for locus_tag in missing_locus_tags]
        logger.error(f"The following genes from the LIMS report format yaml file are missing in the BED file: {', '.join(missing_genes_list)}")
        raise ValueError(f"The following genes from the LIMS report format yaml file are missing in the BED file: {', '.join(missing_genes_list)}")
    else:
        logger.info("All genes from the LIMS report format yaml file are present in the BED file")