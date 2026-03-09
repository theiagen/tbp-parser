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

## Optional Inputs

`tbp-parser` can be customized with a number of optional input parameters. These parameters control:

- files that contain information about the genes of interest and their associated antimicrobials
- files that control the LIMS output report formatting
- quality control thresholds
- text in the output reports (column names, sequencing method, operator name, etc.)

### File Arguments

These options are used to create standard variables that are used throughout the script. These files were previously global variables in earlier versions of tbp-parser, but are now defined as input arguments to allow for greater customization.

| Name | Description | Default Value |
| :------------ | :---------- | :------------ |
| `--config` | the configuration files to use, in YAML format. This argument overrides all other arguments EXCEPT for the other file-type arguments. | |
| `-b`, `--tbdb_bed` | the BED file containing the genes of interest, their locus tags, their associated antimicrobial, and their regions for QC calculations; should be formatted like the TBDB.bed file in TBProfiler | [/data/tbdb.bed](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb.bed) |
| `--lims_report_format_yml` | an optional YAML file that specifies the format of the LIMS report; if not provided, a default format will be used | [/data/default-lims-report-format.yml](https://github.com/theiagen/tbp-parser/blob/main/data/default-lims-report-format.yml) |
| `--gene_database_yml` | an optional YAML file that specifies the gene database information for the genes of interest; if not provided, a default format will be used | [/data/default-gene-database_2026-03-03.yml](https://github.com/theiagen/tbp-parser/blob/main/data/default-gene-database_2026-03-03.yml) |

#### Configuration File

Instead of providing the input parameters on the command line, the ability to provide a configuration file in YAML format is available.

The configuration file will **overwrite all command-line arguments**, except for other file arguments. The configuration file can be provided using the `--config` argument.

To overwrite a variable, please use the find and replace variable in the configuration file. This file is **case-sensitive**.

```yaml
# My laboratory reports "rifampicin" as "rifampin", so I want to rename that text in all of the output files.
# We also use Rv0678 instead of mmpR5 and Rv2983 instead of fbiD, so I want to rename those as well.
# I want an output column renamed in the LIMS report from "Sample Name" to "sample accession"
FIND_AND_REPLACE:
  rifampicin: "rifampin"
  fbiD: "Rv2983"
  mmpR5: "Rv0678"
  "Sample Name": "sample accession"

# I can also overwrite any input parameters, like so. 
# This makes it easy to rerun the same analysis on different 
#  samples without rewriting all of the parameters each time.
MIN_FREQUENCY: 0.1
MIN_PERCENT_LOCI_COVERED: 0.7
TNGS: true
TNGS_FREQUENCY_BOUNDARIES:
- 0.1
- 0.95
TNGS_READ_SUPPORT_BOUNDARIES:
- 100
- 500
```

#### TBDB Bed File

#### LIMS Report Format YAML File

#### Gene Database File

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
| `-fr` | `--find_and_replace` | A JSON string that can be used to specify any text in the output files that should be find-and-replaced with other text. The keys will be the text to find, and the values will be the text to replace it with. This is useful for labs that want to customize the text in their reports (e.g. renaming drugs or genes or output columns). For example, `'{"rifampicin": "rifampin", "fbiD": "Rv2983"}'` | '{}' |

### tNGS-specific Arguments

These options are primarily used for tNGS data.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| `--tngs` | Indicates that the input data was generated using a tNGS protocol. Turns on tNGS-specific features | false |
| `-e`, `--err_bed` | the BED file containing the regions to use for calculating the error rate for tNGS data; should be formatted like the ERR.bed file in TBProfiler | |
| `--use_err_as_brr` | if an ERR BED file is provided, use the ERR regions in place of the TBDB regions for breadth of coverage calculations\nNote: this is an experimental option | |
| `--tngs_frequency_boundaries` | the frequency boundaries (comma-delimited; `lower_f,upper_f`) for tNGS QC reporting, used in conjunction with `--tngs_read_support_boundaries` | 0.1,0.1 |
| `--tngs_read_support_boundaries` | the read support boundaries (comma-delimited; `lower_r,upper_r`) for tNGS QC reporting, used in conjunction with `--tngs_frequency_boundaries` | 10,10 |

Other tNGS arguments are available for backwards compatibility, but they are not recommended for all analyses. See the help message for details.

### Logging Arguments

These options change the verbosity of the `stdout` log.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| `--debug` | The highest level of output verbosity detailing every step of the analysis and logic implemented; overwrites --verbose | false |
