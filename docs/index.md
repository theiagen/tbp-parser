# tbp-parser

!!! warning "Not for Diagnostic Use"
    **CAUTION**: The information produced by this program should **not** be used for clinical reporting unless and until extensive validation has occured in **_your_** laboratory on a stable version. Otherwise, the outputs of tbp-parser are for research use only.

!!! dna "TBProfiler Compatibility"
    **tbp-parser v2+ requires inputs generated with TBProfiler v6.0.0+** Please ensure that you are using the correct version of tbp-parser for your version of TBProfiler.

## Overview

`tbp-parser` is a tool developed in partnership with the California Department of Health (CDPH) to parse the output of [Jody Phelan’s TBProfiler tool](https://github.com/jodyphelan/TBProfiler) into four additional files:

1. A [_Laboratorian_ report](./outputs/laboratorian.md), which contains information about each mutation detected and its associated drug resistance profile in a CSV file.
2. A [_LIMS_ report](./outputs/lims.md), which summarizes the highest severity mutations for each antimicrobial drug and the relevant mutations.
3. A [_Looker_ report](./outputs/looker.md), which condenses the information contained in the Laboratorian report into a format suitable for generating a dashboard in Google’s Looker Studio.
4. A [_coverage_ report](./outputs/coverage.md), which contains the percent breadth of coverage over the minimum depth of each gene relative to the H37Rv reference genome in addition to any warnings, such as any deletions identified in the gene that might have contributed to a reduced percent coverage

Please reach out to us at <support@theiagen.com> if you would like any custom file formats and/or changes to these output files that suit your individual needs.

!!! info "Using tbp-parser as part of TheiaProk?"
    To see the inputs required for `tbp-parser` when run as part of the TheiaProk workflow series, please refer to the documentation for the [TheiaProk](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaprok/) workflow in the Public Health Bioinformatics repository.


## Installation

### Docker

You can use our Docker image to run `tbp-parser` without needing to install any dependencies. To pull the latest version of the Docker image, use the following command:

``` bash
docker pull us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:3.0.0
```

### Locally with Python

`tbp-parser` is availble with `pip`. To install `tbp-parser` and its dependencies, run the following command:

```bash
pip install tbp-parser
```

## Usage

### Example Usage

This shows how the script can be run if used inside the Docker container provided above.

``` text
tbp-parser \
    /path/to/data/tbprofiler_output.json \
    /path/to/data/tbprofiler_output.bam \
    -o "example" \
    --min_depth 12 \
    --min_frequency 0.9 \
    --sequencing_method "Illumina NextSeq" \
    --operator "John Doe"
```

Please note that the BAM file must have the accompanying BAI file in the same directory.
