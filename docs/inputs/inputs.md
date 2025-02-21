---
title: Command-line Arguments
---

The inputs on this page reflect the parameters that are applicable for the command-line tool. To see the inputs required for `tbp-parser` when run as part of the TheiaProk workflow series, please refer to the [TheiaProk Inputs](theiaprok.md) page.

## Required Inputs

`tbp-parser` is designed to run immediately after [Jody Phelan’s TBProfiler tool](https://github.com/jodyphelan/TBProfiler). Only three inputs are required: the JSON file produced by `TBProfiler` and the BAM and BAI file produced by `TBProfiler`.

The JSON file contains information about the mutations detected in the sample: the quality, the type, and if that mutation confers resistance to an antimicrobial drug. The BAM file contains the alignment information for the sample and is needed for determining sequencing quality.

| Parameter  | Description |
| :--------- | :---------- |
| input_json | The path to the JSON file that was produced by `TBProfiler` |
| input_bam  | The path to the BAM file that was produced by `TBProfiler` |

!!! info "BAM index file required"
    The BAM file must have the accompanying BAI file in the same directory. It must also be named exactly the same as the BAM file but ending with a `.bai` suffix.

## Optional Inputs

`tbp-parser` can be customized with a number of optional input parameters. These parameters can be used to control the quality control thresholds, the text that appears in the reports, and the names of the output files. The following is a list of all the input parameters that can be used with `tbp-parser`.

In addition to these arguments, `tbp-parser` also has a `-h, --help` argument that will out the list of possible arguments and their descriptions and a `-v, --version` argument that will print out the version of `tbp-parser` that is installed. Both of these commands exit the program after printing their output.

### Configuration File

Instead of providing the input parameters on the command line, the ability to provide a configuration file in YAML format is available. The configuration file will overwrite all command-line arguments, except for the `--verbose` and `--debug` arguments. The configuration file can be provided using the `--config` argument.

The configuration file can also be used to _overwrite_ the global variables that are in use. The global variables available can be found in the [Global Variables](../algorithm/globals.md) page.

To overwriite a variable, please use the following format in the configuration file. The variable names are case-sensitive.

```yaml
# this variable is found in the globals.py file
GENES_FOR_LIMS:
  - "rpoB"
  - "inhA"
  - "pncA"
  - "inhA"

# although this variable can be set with an input parameter, it must be in uppercase here as it appears in the globals.py file
MIN_DEPTH: 15

# these command-line input parameters are not found in the globals.py file so they are indicated in lowercase
add_cs_lims: True
output_prefix: "Test"
coverage_regions: "/path/to/file"

# only variables that can be found in either globals.py or in the command-line arguments will be used
extra_variable: "This will be ignored"
```

### Quality Control Arguments

These options determine the thresholds for quality control.

| Short Version | Long Version           | Description | Default Value |
| :------------ | :--------------------- | :---------- | :------------ |
| -d            | --min_depth            | The minimum depth of coverage required for a site to pass QC | 10 |
| -c            | --min_percent_coverage | The minimum percentage of a region that has depth above the threshold set by `min_depth` (used for a gene/locus to pass QC) | 100 |
| -s            | --min_read_support     | The minimum read support for a mutation to pass QC | 10 |
| -f            | --min_frequency        | The minimum frequency for a mutation to pass QC (0.1 -> 10%)| 0.1 |
| -r            | --coverage_regions     | A BED file containing the regions to calculate percent coverage for | [/data/tbdb-modified-regions.md](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb-modified-regions.bed) |

### Text Arguments

These options are used verbatim in the reports, or are used to name the output files.

| Short Version | Long Version | Description | Default Value |
| :--- | :--- | :---------- | :------------ |
| -m | --sequencing_method | The sequencing method used to gerneate the data; used in the LIMS & Looker reports. Enclose in quotes if including a space | "Sequencing method not provided" |
| -p | --operator | The operator who ran the analysis; used in the LIMS & Looker reports. Enclose in quotes if including a space | "Operator not provided" |
| -o | --output_prefix | The prefix to use for the output files. Do not include any spaces | "tbp-parser" |

### LIMS Arguments

These options are used to customize the LIMS report

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| --add_cs_lims | Adds Cycloserine (CS) fields to the LIMS report | false |

### tNGS-specific Arguments

These options are primarily used for tNGS data, although all frequency and read support arguments are compatible with WGS data.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| --tngs | Indicates that the input data was generated using the Deeplex + CDPH modified protocol. Turns on tNGS-specific global parameters | false |
| --tngs_expert_regions | A BED file containing the regions to calculate coverage for expert rule regions. This is used to determine coverage quality in the regions where resistance-conferring mutations are found, or where a CDC expert rule is applied. This is not used for QC purposes | [/data/tbdb-expert-regions.bed](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb-expert-regions.bed) |
| --rrs_frequency | The minimum frequency for an _rrs_ mutation to pass QC, as _rrs_ has several problematic sites in the Deeplex tNGS assay | 0.1 |
| --rrl_frequency | The minimum frequency for an _rrl_ mutation to pass QC, as _rrl_ has several problematic sites in the Deeplex tNGS assay | 0.1 |
| --rrs_read_support | The minimum read support for an _rrs_ mutation to pass QC, as _rrs_ has several problematic sites in the Deeplex tNGS assay | 10 |
| --rrl_read_support | The minimum read support for an _rrl_ mutation to pass QC, as _rrl_ has several problematic sites in the Deeplex tNGS assay | 10 |
| --rpob449_frequency | The minimum frequency for an rpoB mutation at protein position 449 to pass QC, as this site is problematic in the Deeplex tNGS assay | 0.1 |
| --etha237_frequency | The minimum frequency for an ethA mutation at protein position 237 to pass QC, as this site is problematic in the Deeplex tNGS assay | 0.1 |

### Logging Arguments

These options change the verbosity of the `stdout` log

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| --verbose | Increases the output verbosity to describe which stage of the analysis is currently running | false |
| --debug | The highest level of output verbosity detailing every step of the analysis and logic implemented; overwrites --verbose | false |
