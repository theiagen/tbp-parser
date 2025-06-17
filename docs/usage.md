---
title: Getting Started
---

## Installation

### Docker

We highly recommend using the following Docker image to run tbp-parser:

``` bash
docker pull us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.6.0 #(1)!
```

1. We host our Docker images on the Google Artifact Registry so that they are always availble for usage.

The entrypoint for this Docker image is the `tbp-parser` help message. To run this container *interactively*, you can use the following command:

``` bash
docker run -it --entrypoint=/bin/bash us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.6.0

# Once inside the container interactively, you can run the tbp-parser tool
python3 /tbp-parser/tbp_parser/tbp_parser.py -v
# 2.6.0
```

### Locally with Python

`tbp-parser` is not yet available with `pip` or `conda`. To run `tbp-parser` in your local command-line environment, install the following dependencies:

- python3
- pandas >= 1.4.2
- importlib_resources
- pyyaml
- samtools

After installation of these dependencies, download and extract the latest release of `tbp-parser` and run the script with `python3`.

## Usage

### Example Usage

This shows how the script can be run if used inside the Docker container provided above.

``` text
python3 /tbp-parser/tbp_parser/tbp_parser.py \
    /path/to/data/tbprofiler_output.json \
    /path/to/data/tbprofiler_output.bam \
    -o "example" \
    --min_depth 12 \
    --min_frequency 0.9 \
    --sequencing_method "Illumina NextSeq" \
    --operator "John Doe" 
```

Please note that the BAM file must have the accompanying BAI file in the same directory.

### Help Message

The help message printed by `tbp-parser` is quite extensive, but has a lot of useful information regarding the input parameters. Here is the entire message in full. You can find more information regarding these inputs in the [Inputs](inputs/inputs.md) section.

``` text
usage: python3 /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]

Parses Jody Phelon's TB-Profiler JSON output into four files:
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
  --config 
          the configuration file to use, in YAML format
          (overrides all other arguments except input_json and input_bam)
          
quality control arguments:
  options that determine what passes QC

  -d, --min_depth
          the minimum depth of coverage for a site to pass QC
          default=10
  -c, --min_percent_coverage
          the minimum percentage of a region that has depth above the threshold set by min_depth
            (used for a gene/locus to pass QC)
          default=100
  -s, --min_read_support
          the minimum read support for a mutation to pass QC
          default=10
  -f, --min_frequency
          the minimum frequency for a mutation to pass QC (0.1 -> 10%)
          default=0.1
  -r, --coverage_regions
          the BED file containing the regions to calculate percent coverage for
          default=data/tbdb-modified-regions.bed

text arguments:
  arguments that are used verbatim in the reports or to name the output files

  -m, --sequencing_method
          the sequencing method used to generate the data; used in the LIMS & Looker reports
          ** Enclose in quotes if includes a space
          default="Sequencing method not provided"
  -p, --operator
          the operator who ran the sequencing; used in the LIMS & Looker reports
          ** Enclose in quotes if includes a space
          default="Operator not provided"
  -o, --output_prefix
          the output file name prefix
          ** Do not include any spaces

tNGS-specific arguments:
  options that are primarily used for tNGS data
  (all frequency arguments are compatible with WGS data)

  --tngs
          indicates that the input data was generated using Deeplex + CDPH modified protocol
          Turns on tNGS-specific global parameters
  --tngs_expert_regions
          the BED file containing the regions to calculate coverage for expert rule regions
            (used to determine coverage quality in the regions where resistance-conferring
            mutations are found, or where a CDC expert rule is applied; not for QC)
          default=data/tngs-expert-rule-regions.bed
  --rrs_frequency
          the minimum frequency for an rrs mutation to pass QC
            (rrs has several problematic sites in the Deeplex tNGS assay)
          default=0.1
  --rrl_frequency
          the minimum frequency for an rrl mutation to pass QC
            (rrl has several problematic sites in the Deeplex tNGS assay)
          default=0.1
  --rpob449_frequency
          the minimum frequency for an rpoB mutation at protein position 449 to pass QC
            (this is a problematic site in the Deeplex tNGS assay)
          default=0.1
  --etha237_frequency
          the minimum frequency for an ethA mutation at protein position 237 to pass QC
            (this is a problematic site in the Deeplex tNGS assay)
          default=0.1

logging arguments:
  options that change the verbosity of the stdout log

  --verbose
          increase output verbosity
  --debug
          increase output verbosity to debug; overwrites --verbose

Please contact support@theiagen.com with any questions
```
