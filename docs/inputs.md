---
title: Command-line Arguments
---

The inputs on this page reflect the parameters that are applicable for the command-line tool. To see the inputs required for `tbp-parser` when run as part of the TheiaProk workflow series, please refer to the documentation for the [TheiaProk](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaprok/) workflow in the Public Health Bioinformatics repository.

## Required Inputs

`tbp-parser` is designed to run immediately after [Jody Phelan’s TBProfiler tool](https://github.com/jodyphelan/TBProfiler). Only three inputs are required: the JSON file produced by `TBProfiler` and the BAM and BAI file produced by `TBProfiler`.

| Parameter  | Description | Purpose |
| :--------- | :---------- | :------ |
| input_json | The path to the JSON file that was produced by `TBProfiler` |Contains information about the mutations detected in the sample: quality, type, and any antimicrobial resistance information. |
| input_bam  | The path to the BAM file that was produced by `TBProfiler` | Contains the alignment information for the sample; needed for determining sequencing quality for quality control. |

!!! info "BAM index file required"
    The BAM file must have the accompanying BAI file in the same directory. It must also be named exactly the same as the BAM file but ending with a `.bai` suffix.


The following global variables are used in the `tbp-parser` algorithm. These variables are used to determine the columns in the LIMS report and to rename certain output text in all output files.

These variables can be overwritten using a configuration file. For more information on how to overwrite these variables, please see the [Configuration File](#configuration-file) section.

## Global Variables

The content of the global variables can be found by browsing the [globals.py](https://github.com/theiagen/tbp-parser/blob/main/tbp_parser/globals.py) file on GitHub. These variables are used to determine the columns in the LIMS report and to rename certain output text in all output files.

These variables can be overwritten using a configuration file. For more information on how to overwrite these variables, please see the [Configuration File](#configuration-file) section.

| Variable Name | Description | Format |
| :------------ | :---------- | :------------ |
| DRUG_COLUMNS_TO_GENE_COLUMNS | This relatively complex dictionary is derived from the [TBProfiler `genes.bed` file](https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb/genes.bed). Used to create the LIMS report, this dictionary can be easily extended to include new drugs/genes or reduced to only the drugs/genes needed in the configuration file. | dict[str, dict[str, dict[str, str]]] |
| OUTPUT_RENAMING | This dictionary maps existing default language to new language for output reports. For example, if the desired output language for 'rifampicin' is 'rifampin', the dictionary would contain the entry {'rifampicin': 'rifampin'}. At the end of processing, all instances of 'rifampicin' in the output reports would be replaced with 'rifampin'. | dict[str, str] |

## Optional Inputs

`tbp-parser` can be customized with a number of optional input parameters. These parameters control:

- quality control thresholds
- miscellaneous text that appears in the reports
- the names of the output files

### Configuration File

Instead of providing the input parameters on the command line, the ability to provide a configuration file in YAML format is available.

The configuration file will **overwrite all command-line or global arguments**, except for the `--verbose` and `--debug` arguments. The configuration file can be provided using the `--config` argument. The global variables available can be found above.

To overwriite a variable, please use the following format in the configuration file. The variable names are **case-sensitive**.

Variables found in the `globals.py` file are written in uppercase, while variables provided as command-line arguments are written in lowercase. This is reflected in the example below.

```yaml
# My laboratory reports "rifampicin" as "rifampin", so I want to rename that text in all of the output files.
OUTPUT_RENAMING:
  rifampicin: "rifampin"

# My laboratory only tests for a subset of drugs, so I want to reduce the output of the LIMS report to only those drugs and genes
# Even though my lab reports rifampicin as rifampin (as indicated in the OUTPUT_RENAMING configuration value above), I need to use the original key here
DRUG_COLUMNS_TO_GENE_COLUMNS:
    rifampicin: 
        MY_UNIQUE_COLUMN_NAME_FOR_RIF:
            rpoB: MY_UNIQUE_COLUMN_NAME_FOR_RPOB
    linezolid:
        MY_UNIQUE_COLUMN_NAME_FOR_LZD:
            rrl: MY_UNIQUE_COLUMN_NAME_FOR_RRL
            rplC: MY_UNIQUE_COLUMN_NAME_FOR_RPLC

# overwrite the following input parameters
min_depth: 15
output_prefix: "Test"
tbdb_bed: "/path/to/file"

# only variables that can be found in either globals.py or in the command-line arguments will be used
extra_variable: "This will be ignored"
IGNORED_FIELD: 12345  # this will also be ignored
```

### File Arguments

These options are used to create standard variables that are used throughout the script. These files were previously global variables in earlier versions of tbp-parser, but are now defined as input arguments to allow for greater customization.

| Short Version | Long Version | Description | Default Value |
| :------------ | :------------ | :---------- | :------------ |
| `-b` | `--tbdb_bed` | the BED file containing the genes of interest, their locus tags, their associated antimicrobial, and their regions for QC calculations; should be formatted like the TBDB.bed file in TBProfiler | [/data/tbdb.bed](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb.bed) |
| `-g` | `--gene_tier_tsv` | a TSV file mapping genes to their tier | [/data/gene-to-tier_2025-12-10.tsv](https://github.com/theiagen/tbp-parser/blob/main/data/gene-to-tier_2025-12-10.tsv) |
| `-p` | `--promoter_regions_tsv` | a TSV file containing the promoter regions to include in interpretation designations; created from the WHO v2 catalogue | [/data/who-v2-promoter-regions_2025-12-10.tsv](https://github.com/theiagen/tbp-parser/blob/main/data/who-v2-promoter-regions_2025-12-10.tsv) |

### Quality Control Arguments

These options determine the thresholds for quality control.

| Short Version | Long Version           | Description | Default Value |
| :------------ | :--------------------- | :---------- | :------------ |
| `-d` | `--min_depth` | The minimum depth of coverage required for a site to pass QC | 10 |
| `-c` | `--min_percent_coverage` | The minimum percentage of a region that has depth above the threshold set by `min_depth` (used for a gene/locus to pass QC; 1.0 -> 100%) | 1.0 |
| `-s` | `--min_read_support` | The minimum read support for a mutation to pass QC | 10 |
| `-f` | `--min_frequency` | The minimum frequency for a mutation to pass QC (0.1 -> 10%) | 0.1 |
| `-l` | `--min_percent_loci_covered` | The minimum percentage of loci/genes in the LIMS report that must pass coverage QC for the sample to be identified as MTBC (0.7 -> 70%) | 0.7 |

### Text Arguments

These options are used verbatim in the reports, or are used to name the output files.

| Short Version | Long Version | Description | Default Value |
| :--- | :--- | :---------- | :------------ |
| `-m` | `--sequencing_method` | The sequencing method used to gerneate the data; used in the LIMS & Looker reports. Enclose in quotes if including a space | "Sequencing method not provided" |
| `-t` | `--operator` | The operator who ran the analysis; used in the LIMS & Looker reports. Enclose in quotes if including a space | "Operator not provided" |
| `-o` | `--output_prefix` | The prefix to use for the output files. Do not include any spaces | "tbp_parser" |

### tNGS-specific Arguments

These options are primarily used for tNGS data.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| `--tngs` | Indicates that the input data was generated using a tNGS protocol. Turns on tNGS-specific features | false |

### Logging Arguments

These options change the verbosity of the `stdout` log.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| `--verbose` | Increases the output verbosity to describe which stage of the analysis is currently running | false |
| `--debug` | The highest level of output verbosity detailing every step of the analysis and logic implemented; overwrites --verbose | false |
