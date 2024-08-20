---
title: Command-line Arguments
---

The inputs on this page reflect the parameters that are applicable for the command-line tool. To see the inputs required for `tbp-parser` when run as part of the TheiaProk workflow series, please refer to the [TheiaProk Inputs](theiaprok.md) page.

## Required Inputs

`tbp-parser` is designed to run immediately after [Jody Phelan’s TB-Profiler tool](https://github.com/jodyphelan/TBProfiler). Only two inputs are required: the JSON file produced by `TB-Profiler` and the BAM file produced by `TB-Profiler`.

The JSON file contains information about the mutations detected in the sample: the quality, the type, and if that mutation confers resistance to an antimicrobial drug. The BAM file contains the alignment information for the sample and is needed for determining sequencing quality. 

| Parameter  | Description |
| :--------- | :---------- |
| input_json | The path to the JSON file that was produced by `TB-Profiler` |
| input_bam  | The path to the BAM file that was produced by `TB-Profiler` |

!!! info
    The BAM file must have the accompanying BAI file in the same directory. It must also be named exactly the same as the BAM file but ending with a `.bai` suffix.

## Optional Inputs

`tbp-parser` can be customized with a number of optional input parameters. These parameters can be used to control the quality control thresholds, the text that appears in the reports, and the names of the output files. The following is a list of all the input parameters that can be used with `tbp-parser`.

In addition to these arguments, `tbp-parser` also has a `-h, --help` argument that will out the list of possible arguments and their descriptions and a `-v, --version` argument that will print out the version of `tbp-parser` that is installed. Both of these commands exit the program after printing their output.

### Quality Control Arguments

These options determine the thresholds for quality control.

| Short Version | Long Version           | Description | Default Value |
| :------------ | :--------------------- | :---------- | :------------ |
| -d            | --min_depth            | The minimum depth of coverage required for a site to pass QC | 10 |
| -c            | --min_percent_coverage | The minimum percentage of a region that has depth above the threshold set by `min_depth` (used for a gene/locus to pass QC) | 100 |
| -s            | --min_read_support     | The minimum read support for a mutation to pass QC | 10
| -f            | --min_frequency        | The minimum frequency for a mutation to pass QC (0.1 -> 10%)| 0.1 |
| -r            | --coverage_regions     | A BED file containing the regions to calculate percent coverage for | [/data/tbdb-modified-regions.md](https://github.com/theiagen/tbp-parser/blob/v1.6.0/data/tbdb-modified-regions.bed) |

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

These options are primarily used for tNGS data, although all frequency arguments are compatible with WGS data.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| --tngs | Indicates that the input data was generated using the Deeplex + CDPH modified protocol. Turns on tNGS-specific global parameters | false |
| --tngs_expert_regions | A BED file containing the regions to calculate coverage for expert rule regions. This is used to determine coverage quality in the regions where resistance-conferring mutations are found, or where a CDC expert rule is applied. This is not used for QC purposes | [/data/tbdb-expert-regions.bed](https://github.com/theiagen/tbp-parser/blob/v1.6.0/data/tbdb-expert-regions.bed) |
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
