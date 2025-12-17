# tbp-parser

!!! warning "Not for Diagnostic Use"
    **CAUTION**: The information produced by this program should **not** be used for clinical reporting unless and until extensive validation has occured in ==_your_== laboratory on a stable version. Otherwise, the outputs of tbp-parser are for research use only.

!!! dna "TBProfiler Compatibility"
    **Versions of TBProfiler prior to v6.0.0 are not compatible with v2+ of tbp-parser.** Please ensure that you are using the correct version of tbp-parser for your version of TBProfiler.

## Overview

`tbp-parser` is a tool developed in partnership with the California Department of Health (CDPH) to parse the output of [Jody Phelan’s TBProfiler tool](https://github.com/jodyphelan/TBProfiler) into four additional files:

1. A _Laboratorian_ report, which contains information about each mutation detected and its associated drug resistance profile in a CSV file.
2. A _LIMS_ report, which summarizes the highest severity mutations for each antimicrobial drug and the relevant mutations.
3. A _Looker_ report, which condenses the information contained in the Laboratorian report into a format suitable for generating a dashboard in Google’s Looker Studio.
4. A _coverage_ report, which contains the percent breadth of coverage over the minimum depth of each gene relative to the H37Rv reference genome in addition to any warnings, such as any deletions identified in the gene that might have contributed to a reduced percent coverage

Please reach out to us at <support@theiagen.com> if you would like any custom file formats and/or changes to these output files that suit your individual needs.

## Installation

### Docker

We highly recommend using the following Docker image to run tbp-parser:

``` bash
docker pull us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:3.0.0
```

The entrypoint for this Docker image is the `tbp-parser` help message. To run this container _interactively_, you can use the following command:

``` bash
docker run -it --entrypoint=/bin/bash us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:3.0.0

# Once inside the container interactively, you can run the tbp-parser tool
python3 /tbp-parser/tbp_parser/tbp_parser.py -v
# 3.0.0
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
