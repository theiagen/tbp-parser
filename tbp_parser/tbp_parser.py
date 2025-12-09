#!/usr/bin/env python3
import argparse
import CheckInputs
import importlib_resources
from __init__ import __VERSION__
from Parser import Parser 

def main():
    home_dir = importlib_resources.files("tbp_parser")
    default_coverage_regions = home_dir.joinpath("..", "data", "tbdb.bed")
    
    parser = argparse.ArgumentParser(
        prog = "tbp-parser",
        description = "Parses Jody Phelon's TBProfiler JSON output into four files:\n- a Laboratorian report,\n- a LIMS report\n- a Looker report, and\n- a coverage report",
        usage = "python3 /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]",
        epilog = "Please contact support@theiagen.com with any questions",
        formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=10))
    parser.add_argument("input_json", 
                        help="the JSON file produced by TBProfiler", type=CheckInputs.is_json_valid)
    parser.add_argument("input_bam", 
                        help="the BAM file produced by TBProfiler", type=CheckInputs.is_bam_valid)
    parser.add_argument("-v", "--version", 
                        action='version', version=str(__VERSION__))
    parser.add_argument("--config",
                        help="the configuration file to use, in YAML format\n(overrides all other arguments except input_json and input_bam)", default="", metavar="\b", type=CheckInputs.is_config_valid)

    qc_arguments = parser.add_argument_group("quality control arguments", 
                                              "options that determine what passes QC")
    qc_arguments.add_argument("-d", "--min_depth", 
                        help="the minimum depth of coverage for a site to pass QC\ndefault=10", default=10, metavar="\b", type=int)
    qc_arguments.add_argument("-c", "--min_percent_coverage", 
                        help="the minimum percentage of a region that has depth above the threshold set by min_depth\n  (used for a gene/locus to pass QC)\ndefault=100", default=100, metavar="\b", type=float)
    qc_arguments.add_argument("-s", "--min_read_support",
                        help="the minimum read support for a mutation to pass QC\ndefault=10", default=10, metavar="\b", type=int)
    qc_arguments.add_argument("-f", "--min_frequency",
                        help="the minimum frequency for a mutation to pass QC (0.1 -> 10%%)\ndefault=0.1", default=0.1, metavar="\b", type=float)
    
    ### TO-DO: consider renaming the following argument to something more intuitive
    qc_arguments.add_argument("-r", "--coverage_regions",
                        help="the BED file containing the regions to calculate percent breadth of coverage for\ndefault=data/tbdb.bed", default=default_coverage_regions, metavar="\b", type=CheckInputs.is_bed_valid)
    
    qc_arguments.add_argument("-l", "--min_percent_locus_covered", default=0.7, metavar="\b", type=float,
                        help="the minimum percentage of loci/genes in the LIMS report that must pass coverage QC for the sample to be identified as MTBC (0.7 -> 70%%)\ndefault=0.7")
    qc_arguments.add_argument("--treat_r_mutations_as_s", default=False, action="store_true",
                        help="treat R mutations the same as S or U mutations in that if locus coverage is poor, they will not be reported regardless of mutation quality\ndefault=False")

    general_arguments = parser.add_argument_group("text arguments", 
                                                  "arguments that are used verbatim in the reports or to name the output files")
    general_arguments.add_argument("-m", "--sequencing_method", 
                        help="the sequencing method used to generate the data; used in the LIMS & Looker reports\n** Enclose in quotes if includes a space\ndefault=\"Sequencing method not provided\"", default="Sequencing method not provided", metavar="\b")
    general_arguments.add_argument("-p", "--operator", 
                        help="the operator who ran the sequencing; used in the LIMS & Looker reports\n** Enclose in quotes if includes a space\ndefault=\"Operator not provided\"", default="Operator not provided", metavar="\b")
    general_arguments.add_argument("-o", "--output_prefix", 
                        help="the output file name prefix\n** Do not include any spaces", default="tbp_parser", metavar="\b")

    tngs_arguments = parser.add_argument_group("tNGS-specific arguments", 
                                                "options that are primarily used for tNGS data\n(all frequency arguments are compatible with WGS data)")
    # TO-DO: reevaluate this argument's help message as it is incorrect
    tngs_arguments.add_argument("--tngs",
                        help="\nindicates that the input data was generated using Deeplex + CDPH modified protocol\nTurns on tNGS-specific global parameters", action="store_true", default=False)

    # TO-DO: reevaluate these arguments -- are they still needed?
    tngs_arguments.add_argument("--rrs_frequency",
                        help="the minimum frequency for an rrs mutation to pass QC\n  (rrs has several problematic sites in the Deeplex tNGS assay)\ndefault=0.1", default=0.1, metavar="\b", type=float)
    tngs_arguments.add_argument("--rrs_read_support",
                        help="the minimum read support for an rrs mutation to pass QC\n  (rrs has several problematic sites in the Deeplex tNGS assay)\ndefault=10", default=10, metavar="\b", type=int)
    tngs_arguments.add_argument("--rrl_frequency",
                        help="the minimum frequency for an rrl mutation to pass QC\n  (rrl has several problematic sites in the Deeplex tNGS assay)\ndefault=0.1", default=0.1, metavar="\b", type=float)
    tngs_arguments.add_argument("--rrl_read_support",
                        help="the minimum read support for an rrl mutation to pass QC\n  (rrl has several problematic sites in the Deeplex tNGS assay)\ndefault=10", default=10, metavar="\b", type=int)
    tngs_arguments.add_argument("--rpob449_frequency",
                        help="the minimum frequency for an rpoB mutation at protein position 449 to pass QC\n  (this is a problematic site in the Deeplex tNGS assay)\ndefault=0.1", default=0.1, metavar="\b", type=float)
    tngs_arguments.add_argument("--etha237_frequency",
                        help="the minimum frequency for an ethA mutation at protein position 237 to pass QC\n  (this is a problematic site in the Deeplex tNGS assay)\ndefault=0.1", default=0.1, metavar="\b", type=float)

    # new arguments for qc reporting
    boundary_arguments = parser.add_argument_group("tNGS-specific QC boundary arguments (NOT compatible with WGS data)",
                                                  "options that set read support and frequency boundaries for tNGS QC reporting:\n if `lower_rs <= read support < upper_rs` AND frequency >= `upper_f` OR\n if `read support >= upper_rs` AND frequency >= `lower_f`\n  QC pass")
    boundary_arguments.add_argument("--tngs_read_support_boundaries",
                        help="the read support boundaries (comma-delimited; \"lower_rs,upper_rs\") for tNGS QC reporting, used in conjunction with `--tngs_frequency_boundaries`\ndefault=10,10 (this is equivalent to no change from the default min read support)",
                        default="10,10", metavar="\b")
    boundary_arguments.add_argument("--tngs_frequency_boundaries",
                        help="the frequency boundaries (comma-delimited; \"lower_f,upper_f\") for tNGS QC reporting, used in conjunction with `--tngs_read_support_boundaries`\ndefault=0.1,0.1 (this is equivalent to no change from the default min frequency)",
                        default="0.1,0.1", metavar="\b")

    logging_arguments = parser.add_argument_group("logging arguments", 
                                                  "options that change the verbosity of the stdout log")
    logging_arguments.add_argument("--verbose", 
                        help="increase output verbosity", action="store_true", default=False)
    logging_arguments.add_argument("--debug", 
                        help="increase output verbosity to debug; overwrites --verbose", action="store_true", default=False)

    options = parser.parse_args()

    parse = Parser(options)
    parse.run()


if __name__ == "__main__":
    main()