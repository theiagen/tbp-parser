#!/usr/bin/env python3
import argparse
import CheckInputs
import importlib_resources
from __init__ import __VERSION__
from Parser import Parser 

def main():
    home_dir = importlib_resources.files("tbp_parser")
    default_coverage_regions = home_dir.joinpath("..", "data", "tbdb.bed")
    default_gene_tier_tsv = home_dir.joinpath("..", "data", "gene-to-tier_2025-12-10.tsv")
    default_promoter_regions = home_dir.joinpath("..", "data", "who-v2-promoters_2025-12-10.tsv")   
    
    parser = argparse.ArgumentParser(
        prog = "tbp-parser",
        description = "Parses Jody Phelon's TBProfiler JSON output into four files:\n- a Laboratorian report,\n- a LIMS report\n- a Looker report, and\n- a coverage report",
        usage = "python3 /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]",
        epilog = "Please contact support@theiagen.com with any questions",
        formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=10))
    parser.add_argument("input_json", 
                        help="the JSON file produced by TBProfiler", type=CheckInputs.is_file_valid)
    parser.add_argument("input_bam", 
                        help="the BAM file produced by TBProfiler", type=CheckInputs.is_bam_valid)
    parser.add_argument("-v", "--version", 
                        action='version', version=str(__VERSION__))
    parser.add_argument("--config",
                        help="the configuration file to use, in YAML format\n(overrides all other arguments EXCEPT for any file-type inputs)", default="", metavar="\b", type=CheckInputs.is_optional_file_valid)

    file_arguments = parser.add_argument_group("file arguments",
                                                "arguments that specify input files used to create standard dictionaries")
    
    ### TO-DO: brainstorm better name for this argument -- regions_bed, targets_bed, genes_bed, etc.
    file_arguments.add_argument("-b", "--tbdb_bed",
                        help="the BED file containing the genes of interest, their locus tags, their associated antimicrobial, and their regions for QC calculations; should be formatted like the TBDB.bed file in TBProfiler\ndefault=data/tbdb.bed", default=default_coverage_regions, metavar="\b", type=CheckInputs.is_bed_valid)
    file_arguments.add_argument("-g", "--gene_tier_tsv", 
                        help="the TSV file mapping genes to their tier\ndefault=data/gene-to-tier_2025-12-10.tsv", default=default_gene_tier_tsv, metavar="\b", type=CheckInputs.is_file_valid)
    file_arguments.add_argument("-p", "--promoter_regions_tsv",
                        help="the TSV file containing the promoter regions to include in interpretation designations\ndefault=data/who-v2-promoters_2025-12-10.tsv", default=default_promoter_regions, metavar="\b", type=CheckInputs.is_file_valid)

    qc_arguments = parser.add_argument_group("quality control arguments", 
                                              "options that determine what passes QC")
    qc_arguments.add_argument("-d", "--min_depth", 
                        help="the minimum depth of coverage for a site to pass QC\ndefault=10", default=10, metavar="\b", type=int)
    qc_arguments.add_argument("-c", "--min_percent_coverage", 
                        help="the minimum percentage of a region that has depth above the threshold set by min_depth\n  (used for a gene/locus to pass QC; 1.0 -> 100%%)\ndefault=1.0", default=1.0, metavar="\b", type=float)
    qc_arguments.add_argument("-s", "--min_read_support",
                        help="the minimum read support for a mutation to pass QC\ndefault=10", default=10, metavar="\b", type=int)
    qc_arguments.add_argument("-f", "--min_frequency",
                        help="the minimum frequency for a mutation to pass QC (0.1 -> 10%%)\ndefault=0.1", default=0.1, metavar="\b", type=float)
    qc_arguments.add_argument("-l", "--min_percent_loci_covered", default=0.7, metavar="\b", type=float,
                        help="the minimum percentage of loci/genes in the LIMS report that must pass coverage QC for the sample to be identified as MTBC (0.7 -> 70%%)\ndefault=0.7")
    
    qc_arguments.add_argument("--do_not_treat_r_mutations_differently", default=False, action="store_true",
                        help="treat R mutations the same as S or U mutations in that if locus coverage is poor, they will not be reported\ndefault=False")

    general_arguments = parser.add_argument_group("text arguments", 
                                                  "arguments that are used verbatim in the reports or to name the output files")
    general_arguments.add_argument("-m", "--sequencing_method", 
                        help="the sequencing method used to generate the data; used in the LIMS & Looker reports\n** Enclose in quotes if includes a space\ndefault=\"Sequencing method not provided\"", default="Sequencing method not provided", metavar="\b")
    general_arguments.add_argument("-t", "--operator", 
                        help="the operator who ran the sequencing; used in the LIMS & Looker reports\n** Enclose in quotes if includes a space\ndefault=\"Operator not provided\"", default="Operator not provided", metavar="\b")
    general_arguments.add_argument("-o", "--output_prefix", 
                        help="the output file name prefix\n** Do not include any spaces", default="tbp_parser", metavar="\b")

    tngs_arguments = parser.add_argument_group("tNGS-specific arguments", 
                                                "options that are primarily used for tNGS data")
    tngs_arguments.add_argument("--tngs",
                        help="\nindicates that the input data was generated using a tNGS protocol\nTurns on tNGS-specific features", action="store_true", default=False)
    tngs_arguments.add_argument("-e", "--err_bed",
                                help="an optional BED file formatted similarly to the --tbdb_bed file but containing ranges that are essential for resistance", default=None, metavar="\b", type=CheckInputs.is_optional_file_valid)

    # new arguments for qc reporting
    boundary_arguments = parser.add_argument_group("tNGS-specific QC boundary arguments (NOT compatible with WGS data)",
                                                  "options that set read support and frequency boundaries for tNGS QC reporting:\n if `lower_rs <= read support < upper_rs` AND frequency >= `upper_f` OR\n if `read support >= upper_rs` AND frequency >= `lower_f`\n  QC pass")
    boundary_arguments.add_argument("--tngs_read_support_boundaries",
                        help="the read support boundaries (comma-delimited; \"lower_rs,upper_rs\") for tNGS QC reporting, used in conjunction with `--tngs_frequency_boundaries`\ndefault=10,10 (this is equivalent to no change from the default min read support)",
                        default="10,10", metavar="\b", type=CheckInputs.is_boundary_valid)
    boundary_arguments.add_argument("--tngs_frequency_boundaries",
                        help="the frequency boundaries (comma-delimited; \"lower_f,upper_f\") for tNGS QC reporting, used in conjunction with `--tngs_read_support_boundaries`\ndefault=0.1,0.1 (this is equivalent to no change from the default min frequency)",
                        default="0.1,0.1", metavar="\b", type=CheckInputs.is_boundary_valid)

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