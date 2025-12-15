import os
import argparse

def is_file_valid(filename) -> str:
    """Checks if an input file is accessible

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    if not os.path.exists(filename) and filename != "-":
        raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    return filename

def is_optional_file_valid(filename) -> str:
    """Checks if an optional input file is accessible (no default file provided)

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    if filename != "":
        if not os.path.exists(filename) and filename != "-":
            raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    return filename

def is_bam_valid(filename) -> str:
    """Checks if the input BAM is valid and accessible and if there is an associated BAI

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    if not os.path.exists(filename) and (filename != "-" and not filename.endswith(".bam")):
        raise argparse.ArgumentTypeError("{0} cannot be accessed or is missing the BAM extension".format(filename))
    bai_filename = filename + ".bai"
    if not os.path.exists(bai_filename) and (bai_filename != "-"):
        raise argparse.ArgumentTypeError("Cannot find the BAI file for this BAM")
    return filename

def is_bed_valid(filename) -> str:
    """Checks if the TBDB BED file is accessible

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """

    # check if the necessary columns are present in the BED file -- just count them because we can't really parse it here
    # does this file have at least 6 columns
    if filename != "" and not os.path.exists(filename) and filename != "-":
        raise argparse.ArgumentTypeError("{0} cannot be accessed".format(filename))
    else:
        with open(filename, 'r') as bed_file:
            for line in bed_file:
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    raise argparse.ArgumentTypeError("{0} does not have at least 6 columns as required".format(filename))
                break  # only need to check the first line
    return filename
