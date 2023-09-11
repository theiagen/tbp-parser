#!/usr/bin/env python3
import argparse
import CheckInputs
from __init__ import __VERSION__
from Parser import Parser 

def main():
  parser = argparse.ArgumentParser(
    description = "Parses Jody Phelon's TBProfiler JSON output into three files:\n- a Laboratorian report,\n- a LIMS report\n- a Looker report, and\n- a coverage report",
    usage = "tbp_parser [-h|-v] <input_json> <input_bam> [<args>]",
    formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=10))
  parser.add_argument("input_json", 
                      help="the JSON file produced by TBProfiler", type=CheckInputs.is_json_valid)
  parser.add_argument("input_bam", 
                      help="the BAM file produced by TBProfiler", type=CheckInputs.is_bam_valid)
  parser.add_argument("-v", "--version", 
                      action='version', version=str(__VERSION__))
  parser.add_argument("-o", "--output_prefix", 
                      help="the output file name prefix\nDo not include a space", default="tbprofiler_parser", metavar="\b")
  parser.add_argument("-d", "--min_depth", 
                      help="the minimum depth of coverage to pass QC\ndefault=10", default=10, metavar="\b")
  parser.add_argument("-c", "--coverage_threshold", 
                      help="the minimum percent coverage for a gene to pass QC\ndefault=100", default=100, metavar="\b")
  parser.add_argument("-s", "--sequencing_method", "-s", 
                      help="the sequencing method used to generate the data\nEnclose in quotes if includes a space\ndefault=\"Sequencing method not provided\"", default="Sequencing method not provided", metavar="\b")
  parser.add_argument("-p", "--operator", 
                      help="the operator who ran the sequencing\nEnclose in quotes if includes a space\ndefault=\"Operator not provided\"", default="Operator not provided", metavar="\b")
  parser.add_argument("--verbose", 
                      help="increase output verbosity", action="store_true", default=False)
  parser.add_argument("--debug", 
                      help="increase output verbosity to debug; overwrites --verbose", action="store_true", default=False)

  options = parser.parse_args()
  
  parse = Parser(options)
  parse.run()
  
  
if __name__ == "__main__":
  main()