# tbprofiler_parser

This repository contains the tbprofiler-parser tool which parses the JSON output of [Jody Phelan's TBProfiler tool](https://github.com/jodyphelan/TBProfiler). Available as a download-able Python package and as a Docker image, tbprofiler-parser converts the output of TBProfiler into four files:

1. A _Laboratorian_ report, which contains information regarding each mutation and its associated drug resistance profile in a CSV file. This file also contains two interpretation fields -- "Looker" and "MDL" which are generated using the CDC's expert rules for interpreting the severity of potential drug resistance mutations.
2. A _LIMS_ report, formatted specifically for STAR LIMS. This CSV report summarizes the highest severity mutations for each antimicrobial and lists the relevant mutations for each gene.
3. A _Looker_ report that contains condensed information suitable for generating a dashboard in Google's Looker studio.
4. A _coverage_ report that contains the percent coverage of each gene relative to the H37Rv reference genome in addition to any warnings, such as any deletions identified in the gene.

Please reach out to us at [theiagen@support.com](mailto:theiagen@support.com) if you would like any custom file formats and/or changes to these output files that suit your individual needs.

## Installation

### Docker

We recommend using the following Docker image to run tbprofiler-parser:

To-do: make a docker image

```bash
docker pull theiagen/tbprofiler-parser
```

### Pip/Conda

To-do: add instructions

## Usage

```bash
usage: tbprofiler_parser [-h|-v] <input_json> <input_bam> [<args>]

Parses Jody Phelon's TBProfiler JSON output into three files:
- a Laboratorian report,
- a LIMS report
- a Looker report, and
- a coverage report

positional arguments:
  input_json
          the JSON file produced by TBProfiler
  input_bam
          the BAM file produced by TBProfiler

optional arguments:
  -h, --help
          show this help message and exit
  -v, --version
          show program's version number and exit
  -o, --output_prefix 
          the output file name prefix
          Do not include a space
  -d, --min_depth 
          the minimum depth of coverage to pass QC
          default=10
  -c, --coverage_threshold 
          the minimum percent coverage for a gene to pass QC
          default=100
  -s, --sequencing_method, -s 
          the sequencing method used to generate the data
          Enclose in quotes if includes a space
          default="Sequencing method not provided"
  -p, --operator 
          the operator who ran the sequencing
          Enclose in quotes if includes a space
          default="Operator not provided"
  --verbose
          increase output verbosity
  --debug
          increase output verbosity to debug;
          overwrites --verbose
```
