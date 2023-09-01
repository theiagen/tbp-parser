#!/usr/bin/env python3
import sys
import argparse
import os
import pkg_resources
from tbprofiler_parser.InputTypes import InputTypes
from tbprofiler_parser.Parser import Parser

VERSION = "0.0.1"

parser = argparse.ArgumentParser(
  description='Parses the TBProfiler output into three files for CDPH',
  usage='tbprofiler_parser <input_json> <input_bam> [<args>]'
)
parser.add_argument("input_json", help="the JSON file produced by TBProfiler", type=InputTypes.is_json_valid)
parser.add_argument("input_bam", help="the BAM file produced by TBProfiler", type=InputTypes.is_bam_valid)
parser.add_argument("--verbose", "-v", help="increase output verbosity", action="store_true", default=False)
parser.add_argument("--version", action='version', version=str(VERSION))
parser.add_argument("--output_prefix", "-o", help="the output file name prefix", default="tbprofiler_parser")
parser.add_argument("--min_depth", "-d", help="the minimum depth of coverage to pass QC", default=10)
parser.add_argument('--coverage_threshold, -c', help="the minimum percent coverage for a gene to pass QC", default=100)

options = parser.parse_args()


def tbprofiler_parse(input_json, input_bam):
  """The main script for the TBProfiler_Parser python tool.
  Requires as input the output from Jody Phelan's TBProfiler tool.
  
  Parameters
  ----------
  input_json : `file`
    the JSON file produced by TBProfiler
  input_bam : `file`
    the BAM file produced by TBProfiler
  
  Returns
  -------
  laboratorian_report : `file`
    the CDPH Laboratorian Report
  looker_report : `file`
    the CDPH Looker Report
  lims_report : `file`
    the CDPH LIMS Report
  
  Notes
  -----
  Created for CDPH. Customization available upon request.
  """
  parse = Parser(options)
  parse.run()
  
  
  