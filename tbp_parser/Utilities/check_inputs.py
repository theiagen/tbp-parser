import argparse
import os
import logging
import pysam
from pathlib import Path

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

def _touch_stale_index(filename: str) -> None:
    """
    Bumps a BAM index's mtime if it's older than the BAM itself.

    Removes an unnecessary htslib warning that the index file is older than the data file.
    Touching the index file is computationally simpler than regenerating it.
    """
    bam_path = Path(filename)
    for file in (Path(str(bam_path) + ".bai"), bam_path.with_suffix(".bai")):
        if file.exists():
            if file.stat().st_mtime < bam_path.stat().st_mtime:
                file.touch()
            return

def is_bam_index_valid(filename: str) -> str:
    """Checks if there's an associated BAI for the BAM

    Args:
        filename (String): The name of file to check

    Returns:
        String: The name of the file if valid and accessible
    """
    _touch_stale_index(filename)

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

def validate_err_coords(bed_records, err_records) -> None:
    """
    Checks that every ERR region falls within the target region it belongs to.

    `TargetCoverage.check_err_within_coords` enforces the same invariant, but raises on the first bad
    region and only after the BAM pileup has run. Checking here reports every offending region at
    once, before any BAM work is done.

    Targets are keyed by gene name to mirror `generate_coverage_maps`, which builds one
    TargetCoverage per gene name from that record's coords.

    Args:
        bed_records: List of BedRecord objects from the `--coverage_bed` file
        err_records: List of BedRecord objects from the `--err_coverage_bed` file

    Raises:
        ValueError: If any ERR region is not contained by the target region of the same gene
    """
    target_coords = {record.gene_name: record.coords for record in bed_records}

    outside = []
    for err in err_records:
        coords = target_coords.get(err.gene_name)
        # an ERR region with no matching target is ignored by the coverage calculator, not an error
        if coords is None:
            continue

        start, end = coords
        if not (start <= err.start and err.end <= end):
            outside.append(f"{err.gene_name} ({err.locus_tag}) {err.coords} falls outside {coords}")

    if outside:
        message = "The following ERR regions fall outside their target regions:\n  " + "\n  ".join(outside)
        logger.error(message)
        raise ValueError(message)

def _check_bed_against_gene_db(bed_records) -> list:
    """
    Checks that every locus tag in the BED file exists in the gene database.
    NOTE: The BED file's drug column is parsed into BedRecord.drugs but deliberately not checked
    here. Drug column only parsed from the `--db_bed` file in `build_gene_db`. The truth set of gene/drug interactions
    should already be established at this point in the GeneDatabase.

    Args:
        bed_records: List of BedRecord objects to check

    Returns:
        list: Error messages, empty when every locus tag resolves
    """
    unknown = sorted(
        set(rec.locus_tag for rec in bed_records if GeneDatabase.get_locus_tag(rec.locus_tag) is None)
    )
    if not unknown:
        return []
    return [f"The following genes from the BED file are missing in the Gene Database: {', '.join(unknown)}"]

def _check_lims_against_gene_db(lims_records) -> list:
    """
    Checks the LIMS report format yml file against the gene database for the following error conditions:
      - a drug no gene in the gene database is associated with.
      - a gene code the gene database cannot resolve, whether written as a gene name, locus tag, or alias.
      - a gene/drug pairing the gene database does not contain.

    Args:
        lims_records: List of LIMSRecord objects to check

    Returns:
        list: Error messages, empty when the LIMS entries are a subset of the gene database
    """
    known_drugs = {drug for entry in GeneDatabase.get_db().values() for drug in (entry.get("drugs") or [])}

    unknown_drugs = set()
    unknown_genes = set()
    missing_pairs = set()

    for rec in lims_records:
        drug_is_known = rec.drug in known_drugs
        if not drug_is_known:
            unknown_drugs.add(rec.drug)

        for gene in rec.gene_codes.keys():
            locus_tag = GeneDatabase.get_locus_tag(gene)
            if locus_tag is None:
                unknown_genes.add(gene)
                continue

            # an unknown drug is already reported on its own
            if drug_is_known and rec.drug not in GeneDatabase.get_drugs(locus_tag):
                missing_pairs.add(f"{rec.drug}|{gene}|{locus_tag}")

    errors = []
    if unknown_drugs:
        errors.append(f"The following drugs from the LIMS report format yaml file are missing in the Gene Database: {', '.join(sorted(unknown_drugs))}")
    if unknown_genes:
        errors.append(f"The following genes from the LIMS report format yaml file are missing in the Gene Database: {', '.join(sorted(unknown_genes))}")
    if missing_pairs:
        errors.append(f"The following drug|gene|locus_tag associations from the LIMS report format yaml file are missing in the Gene Database: {', '.join(sorted(missing_pairs))}")

    return errors

def _check_lims_against_bed(lims_records, bed_records) -> list:
    """
    Checks that every gene in the LIMS report format yml file has a region in the BED file.
    Without a BED region there is no coverage figure for the gene and the LIMS coverage lookup raises.

    Args:
        lims_records: List of LIMSRecord objects to check
        bed_records: List of BedRecord objects the LIMS genes must have a region in

    Returns:
        list: Error messages, empty when every LIMS gene has a BED region
    """
    bed_locus_tags = {rec.locus_tag for rec in bed_records}

    missing = set()
    unresolved = set()

    for rec in lims_records:
        for gene in rec.gene_codes.keys():
            locus_tag = GeneDatabase.get_locus_tag(gene)
            if locus_tag is None:
                unresolved.add(gene)
            elif locus_tag not in bed_locus_tags:
                missing.add(f"{gene}|{locus_tag}")

    errors = []
    if missing:
        errors.append(f"The following gene|locus_tag entries from the LIMS report format yaml file are missing in the BED file: {', '.join(sorted(missing))}")
    if unresolved:
        errors.append(f"The following genes from the LIMS report format yaml file could not be checked against the BED file because they do not resolve in the Gene Database: {', '.join(sorted(unresolved))}")

    return errors

def _check_variants_against_gene_db(variant_records) -> list:
    """
    Checks that every gene/drug pairing in the input results JSON exists in the gene database.

    The gene database must be a faithful representation of the TBProfiler database the variants were
    called against, so the results JSON can only ever be a subset of it. A gene or drug the JSON
    reports on but the gene database does not know is potentially reporting misleading information downstream
    because of `_expand_annotations_for_all_drugs` and `generate_unreported_variants`.

    Args:
        variant_records: List of VariantRecord objects to check

    Returns:
        list: Error messages, empty when every gene and pairing is known
    """
    unknown_genes = set()
    missing_pairs = set()

    for rec in variant_records:
        locus_tag = GeneDatabase.get_locus_tag(rec.gene_id)
        if locus_tag is None:
            unknown_genes.add(rec.gene_id)
            continue

        # truth set of drugs from the gene_db
        known_drugs = GeneDatabase.get_drugs(locus_tag)

        # drugs reported via the input results JSON (annotation + gene_associated_drugs)
        annotated_drugs = {annotation.drug for annotation in rec.annotation if annotation.drug}
        all_drugs = annotated_drugs | set(rec.gene_associated_drugs)

        for drug in all_drugs:
            if drug not in known_drugs:
                missing_pairs.add(f"{drug}|{rec.gene_id}|{locus_tag}")

    errors = []
    if unknown_genes:
        errors.append(f"The following genes from the results JSON file are missing in the Gene Database: {', '.join(sorted(unknown_genes))}")
    if missing_pairs:
        errors.append(f"The following drug|gene|locus_tag associations from the results JSON file are missing in the Gene Database: {', '.join(sorted(missing_pairs))}")

    return errors

def validate_inputs(
    bed_records,
    lims_records,
    variant_records,
) -> None:
    """
    Checks that the BED, LIMS report format, and results JSON inputs are all subsets of the gene database.

    The gene database represents the upstream TBProfiler database that variants were called against,
    so anything in the other input files that is missing from it can never be reported on. Input
    files that cover fewer genes or drugs than the gene database are fine; extras are not.

    Every check runs before anything is raised, so a broken set of input files can be corrected in a
    single pass rather than one error at a time.

    Args:
        bed_records: List of BedRecord objects to check
        lims_records: List of LIMSRecord objects to check
        variant_records: Optional list of VariantRecord objects from the results JSON file

    Raises:
        ValueError: If any input file references a gene, drug, or gene/drug pairing that the gene
                    database does not contain, or a LIMS gene with no region in the BED file
    """
    errors = [
        *_check_variants_against_gene_db(variant_records),
        *_check_bed_against_gene_db(bed_records),
        *_check_lims_against_gene_db(lims_records),
        *_check_lims_against_bed(lims_records, bed_records),
    ]

    if errors:
        message = "\n".join(errors) + "\nEither correct the input files, or provide a --gene_database_yml file that contains these entries (see `tbp-parser build_gene_db`)."
        logger.error(message)
        raise ValueError(message)

    logger.info("All genes and gene/drug associations from the input files are present in the Gene Database.")