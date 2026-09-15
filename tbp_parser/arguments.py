#!/usr/bin/env python3
import argparse
import json
import logging
import re
from pathlib import Path
from tbp_parser import __VERSION__
from tbp_parser.Utilities import (
  is_boundary_valid,
  is_fraction_valid,
  is_file_valid,
  is_bam_index_valid,
  is_optional_file_valid,
)

logger = logging.getLogger(__name__)

BUILD_GENE_DB_COMMAND = "build_gene_db"
BUILD_LIMS_FMT_COMMAND = "build_lims_fmt"
PARSE_COMMAND = "parse"

def resolve_output_prefix(output_prefix: str) -> Path:
    """
    Resolve the output prefix to handle directories and file names.
    Give default file name if only a directory is provided.
    """
    output_path = Path(output_prefix)
    if output_prefix.endswith('/'):
        output_path.mkdir(parents=True, exist_ok=True)
        output_path = output_path / "tbp_parser"
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path

class CustomFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Custom help formatter that:
        - preserves manual '\n' line breaks in parser/group descriptions (inherited)
        - never shows a metavar/placeholder after the flag name
        - widens the flag column so help text starts at a consistent position
    """

    def __init__(self, prog, indent_increment=2, max_help_position=40):
        super().__init__(prog, indent_increment, max_help_position)

    def _get_default_metavar_for_optional(self, action):
        # used to build the top usage line cleanly
        return ''

    def _format_actions_usage(self, actions, groups):
        # fixes the trailing space left from removing the metavar in the usage section
        text = super()._format_actions_usage(actions, groups)
        return re.sub(r" {2,}", " ", text).strip()

def parse_arguments(argv: list | None = None):

    main = argparse.ArgumentParser(
        prog = "tbp-parser",
        description = "A TBProfiler output parser used to classify the severity of Mycobacterium tuberculosis mutations.",
        epilog = "Please contact support@theiagen.com with any questions",
        formatter_class = CustomFormatter,
    )
    main.add_argument("-v", "--version", action='version', version=str(__VERSION__))
    subcommands = main.add_subparsers(title="SUBCOMMANDS", dest="command", required=True)

    # =========================================================================
    # PARSE - SUBCOMMAND
    # =========================================================================
    parser = subcommands.add_parser(
        PARSE_COMMAND,
        help="Parse a TBProfiler results JSON into the Laboratorian, LIMS, Looker, and coverage reports.",
        description="Parses Jody Phelon's TBProfiler JSON output into four files:\n- a Laboratorian report,\n- a LIMS report\n- a Looker report, and\n- a coverage report",
        formatter_class = CustomFormatter,
    )

    required_files = parser.add_argument_group("REQUIRED")
    required_files.add_argument("--input_json", help="The results JSON file produced by TBProfiler.", required=True, type=is_file_valid)
    required_files.add_argument("--input_bam", help="The BAM file produced by TBProfiler.", required=True, type=is_bam_index_valid)
    required_files.add_argument("--coverage_bed", help="A BED file containing genes of interest and their coordinates.", required=True, type=is_file_valid)
    required_files.add_argument("--lims_report_format_yml", help=f"A YAML file that defines the formatting and structure of the final LIMS report output. See subcommand: `{BUILD_LIMS_FMT_COMMAND}`.", required=True, type=is_file_valid)
    required_files.add_argument("--gene_database_yml", help=f"A YAML file that defines a gene database and represents the TBDB database used to generate the results JSON. See subcommand: `{BUILD_GENE_DB_COMMAND}`.", required=True, type=is_file_valid)

    validation = parser.add_argument_group("VALIDATION")
    validation.add_argument("--skip_input_validation", help="Skip validation that all genes and gene/drug pairs referenced in the input files are present in the gene database.", action="store_true", default=False)

    qc_arguments = parser.add_argument_group("THRESHOLDS")
    qc_arguments.add_argument("--min_depth", help="The minimum depth of coverage for a position to pass QC (default: %(default)s).", default=10, type=int)
    qc_arguments.add_argument("--min_percent_coverage", help="The minimum fraction of a loci/gene that must meet `--min_depth` for that loci/gene to pass coverage QC (default: %(default)s).", default=1.0, type=is_fraction_valid)
    qc_arguments.add_argument("--min_read_support", help="The minimum number of reads supporting a mutation for it to pass QC (default: %(default)s).", default=10, type=int)
    qc_arguments.add_argument("--min_frequency", help="The minimum frequency for a mutation to pass QC (default: %(default)s).", default=0.1, type=is_fraction_valid)
    qc_arguments.add_argument("--min_percent_loci_covered", help="The minimum fraction of loci/genes in the LIMS report that must pass coverage QC for the sample to be identified as MTBC (default: %(default)s).", default=0.7, type=is_fraction_valid)

    tngs_arguments = parser.add_argument_group("tNGS")
    tngs_arguments.add_argument("--tngs", help="Indicate that the input data was generated with a tNGS protocol; turns on tNGS-specific features.", action="store_true", default=False)
    tngs_arguments.add_argument("--err_coverage_bed", help="An optional BED file, formatted like `--coverage_bed`, containing the regions that are essential for resistance.", default=None, type=is_optional_file_valid)
    tngs_arguments.add_argument("--use_err_for_qc", help="Use the `--err_coverage_bed` regions in place of the `--coverage_bed` regions for breadth of coverage calculations. Experimental.", action="store_true", default=False)
    tngs_arguments.add_argument("--resolve_overlapping_regions", help="Resolve overlapping BED regions so that reads shared by overlapping targets are not counted twice. Recommended for tNGS amplicon data.", action="store_true", default=False)

    boundary_arguments = parser.add_argument_group("tNGS THRESHOLDS", "Options that provide additional criteria under which a mutation may pass QC beyond the standard --min_read_support and --min_frequency thresholds.\nA mutation passes QC if (lower_rs <= read support < upper_rs) AND (frequency >= upper_f), or if (read support >= upper_rs) AND (frequency >= lower_f).")
    boundary_arguments.add_argument("--tngs_read_support_boundaries", help="The read support boundaries for tNGS QC reporting, comma-delimited as \"lower_rs,upper_rs\". Used together with `--tngs_frequency_boundaries`. The default is equivalent to applying `--min_read_support` alone (default: %(default)s).", default="10,10", type=is_boundary_valid)
    boundary_arguments.add_argument("--tngs_frequency_boundaries", help="The frequency boundaries for tNGS QC reporting, comma-delimited as \"lower_f,upper_f\". Used together with `--tngs_read_support_boundaries`. The default is equivalent to applying `--min_frequency` alone (default: %(default)s).", default="0.1,0.1", type=is_boundary_valid)

    general_arguments = parser.add_argument_group("TEXT")
    general_arguments.add_argument("--sequencing_method", help="The sequencing method used to generate the data; written verbatim to the LIMS and Looker reports. Enclose in quotes if it contains a space.", default="Sequencing method not provided")
    general_arguments.add_argument("--operator", help="The operator who ran the sequencing; written verbatim to the LIMS and Looker reports. Enclose in quotes if it contains a space.", default="Operator not provided")
    general_arguments.add_argument("--output_prefix", help="The prefix for all output file names; a trailing slash is treated as a directory. Do not include spaces (default: %(default)s).", default="tbp_parser", type=resolve_output_prefix)
    general_arguments.add_argument("--find_and_replace", help="A JSON string of find-and-replace pairs applied across all fields. Example: --find_and_replace '{\"rifampicin\": \"rifampin\", \"fbiD\": \"Rv2983\"}'", default={}, type=json.loads)

    supplementary = parser.add_argument_group("SUPPLEMENTARY")
    supplementary.add_argument("--config", help="The YAML configuration file to use. Values in this file override all other arguments, except file-type inputs.", type=is_optional_file_valid)

    logging_arguments = parser.add_argument_group("LOGGING")
    logging_arguments.add_argument("--debug", help="Increase output verbosity to debug.", action="store_true", default=False)

    # =========================================================================
    # BUILD_GENE_DB - SUBCOMMAND
    # =========================================================================
    build_gene_db = subcommands.add_parser(
        BUILD_GENE_DB_COMMAND,
        help="Build a `--gene_database_yml` file from the 'genes.bed' file describing the TBProfiler database its variants were called against.",
        formatter_class = CustomFormatter,
    )
    build_gene_db.add_argument("--db_bed", help="A BED file containing every gene TBProfiler calls variants against, its locus tag, and the drugs it is associated with.", required=True, type=is_file_valid)
    build_gene_db.add_argument("--output", help="The path to write the gene database YAML file to.", required=True)
    build_gene_db.add_argument("--debug", help="Increase output verbosity to debug.", action="store_true", default=False)

    # =========================================================================
    # BUILD_LIMS_FMT - SUBCOMMAND
    # =========================================================================
    build_lims_fmt = subcommands.add_parser(
      BUILD_LIMS_FMT_COMMAND,
      help="Build a default `--lims_report_format_yml` file based on the gene-drug interactions found in a `--gene_database_yml` file.",
      formatter_class = CustomFormatter,
    )

    build_lims_fmt.add_argument("--gene_database_yml", help=f"The gene database YAML file to derive the LIMS report format from; either user-provided or generated by `tbp-parser {BUILD_GENE_DB_COMMAND}`.", required=True, type=is_file_valid)
    build_lims_fmt.add_argument("--output", help="The path to write the LIMS report format YAML file to.", required=True)
    build_lims_fmt.add_argument("--debug", help="Increase output verbosity to debug.", action="store_true", default=False)

    # collect all arguments
    options = main.parse_args(argv)
    return options