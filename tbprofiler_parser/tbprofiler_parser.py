#!/usr/bin/env python3
import argparse
import CheckInputs
from Parser import Parser 

VERSION = "0.0.1"

parser = argparse.ArgumentParser(
  description='Parses the TBProfiler output into three files for CDPH',
  usage='tbprofiler_parser <input_json> <input_bam> [<args>]')
parser.add_argument("input_json", help="the JSON file produced by TBProfiler", type=CheckInputs.is_json_valid)
parser.add_argument("input_bam", help="the BAM file produced by TBProfiler", type=CheckInputs.is_bam_valid)
parser.add_argument("--verbose", "-v", help="increase output verbosity", action="store_true", default=False)
parser.add_argument("--version", action='version', version=str(VERSION))
parser.add_argument("--output_prefix", "-o", help="the output file name prefix", default="tbprofiler_parser")
parser.add_argument("--min_depth", "-d", help="the minimum depth of coverage to pass QC", default=10)
parser.add_argument("--coverage_threshold", "-c", help="the minimum percent coverage for a gene to pass QC", default=100)
parser.add_argument("--sequencing_method", "-s", help="the sequencing method used to generate the data", default="Sequencing method not provided")
parser.add_argument("--operator", "-p", help="the operator who ran the sequencing", default="Operator not provided")

options = parser.parse_args()

parse = Parser(options)
parse.run()