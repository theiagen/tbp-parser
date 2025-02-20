#!/usr/bin/env python3
import argparse
import CheckInputs
from __init__ import __VERSION__
from Parser import Parser 

def main():
  parser = argparse.ArgumentParser(
    prog = "tbp-parser",
    description = "Parses Jody Phelon's TBProfiler JSON output into three files:\n- a Laboratorian report,\n- a LIMS report\n- a Looker report, and\n- a coverage report",
    usage = "python3 /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]",
    epilog = "Please contact support@theiagen.com with any questions",
    formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=10))
  parser.add_argument("input_json", 
                      help="the JSON file produced by TBProfiler", type=CheckInputs.is_json_valid)
  parser.add_argument("input_bam", 
                      help="the BAM file produced by TBProfiler", type=CheckInputs.is_bam_valid)
  parser.add_argument("-v", "--version", 
                      action='version', version=str(__VERSION__))
  
  qc_arguments = parser.add_argument_group("quality control arguments", "options that determine what passes QC")
  qc_arguments.add_argument("-d", "--min_depth", 
                      help="the minimum depth of coverage for a site to pass QC\ndefault=10", default=10, metavar="\b", type=int)
  qc_arguments.add_argument("-c", "--min_percent_coverage", 
                      help="the minimum percentage of a region that has depth above the threshold set by min_depth\n  (used for a gene/locus to pass QC)\ndefault=100", default=100, metavar="\b", type=int)
  qc_arguments.add_argument("-s", "--min_read_support",
                      help="the minimum read support for a mutation to pass QC\ndefault=10", default=10, metavar="\b", type=int)
  qc_arguments.add_argument("-f", "--min_frequency",
                      help="the minimum frequency for a mutation to pass QC (0.1 -> 10%%)\ndefault=0.1", default=0.1, metavar="\b", type=float)
  qc_arguments.add_argument("-r", "--coverage_regions",
                      help="the BED file containing the regions to calculate percent coverage for\ndefault=data/tbdb-modified-regions.bed", default="../data/tbdb-modified-regions.bed", metavar="\b", type=CheckInputs.is_bed_valid)
  
  general_arguments = parser.add_argument_group("text arguments", "arguments that are used verbatim in the reports or to name the output files")
  general_arguments.add_argument("-m", "--sequencing_method", 
                      help="the sequencing method used to generate the data; used in the LIMS & Looker reports\n** Enclose in quotes if includes a space\ndefault=\"Sequencing method not provided\"", default="Sequencing method not provided", metavar="\b")
  general_arguments.add_argument("-p", "--operator", 
                      help="the operator who ran the sequencing; used in the LIMS & Looker reports\n** Enclose in quotes if includes a space\ndefault=\"Operator not provided\"", default="Operator not provided", metavar="\b")
  general_arguments.add_argument("-o", "--output_prefix", 
                      help="the output file name prefix\n** Do not include any spaces", default="tbp_parser", metavar="\b")
  
  lims_arguments = parser.add_argument_group("LIMS arguments", "options that are used to customize the LIMS report")
  lims_arguments.add_argument("--add_cs_lims",
                              help="adds cycloserine (CS) fields to the LIMS report", action="store_true", default=False)
  
  tngs_arguments = parser.add_argument_group("tNGS-specific arguments", "options that are primarily used for tNGS data\n(all frequency arguments are compatible with WGS data)")
  tngs_arguments.add_argument("--tngs",
                      help="\nindicates that the input data was generated using Deeplex + CDPH modified protocol\nTurns on tNGS-specific global parameters", action="store_true", default=False)
  tngs_arguments.add_argument("--tngs_expert_regions",
                      help="the BED file containing the regions to calculate coverage for expert rule regions\n  (used to determine coverage quality in the regions where resistance-conferring\n  mutations are found, or where a CDC expert rule is applied; not for QC)\ndefault=data/tngs-expert-rule-regions.bed", default="", metavar="\b", type=CheckInputs.is_bed_valid) #../data/tngs-expert-rule-regions.bed
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
  
  logging_arguments = parser.add_argument_group("logging arguments", "options that change the verbosity of the stdout log")
  logging_arguments.add_argument("--verbose", 
                      help="increase output verbosity", action="store_true", default=False)
  logging_arguments.add_argument("--debug", 
                      help="increase output verbosity to debug; overwrites --verbose", action="store_true", default=False)

  options = parser.parse_args()
  
  parse = Parser(options)
  parse.run()
  
  
if __name__ == "__main__":
  main()