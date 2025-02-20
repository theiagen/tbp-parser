---
title: TheiaProk Inputs on Terra
---

When running `tbp-parser` as part of the TheiaProk workflow series ([find documentation for TheiaProk here](https://theiagen.github.io/public_health_bioinformatics/latest/workflows/genomic_characterization/theiaprok/)) on [Terra.bio](https://terra.bio), an optional input must be activated to instruct TheiaProk to run `tbp-parser`.

`tbp-parser` is not on by default due to the nature of this tool and its outputs.

!!! info annotate "TheiaProk Version"
    This information only corresponds to the upcoming PHB v2.3.0 release. These inputs and outputs may not be applicable to other versions of TheiaProk.

*[PHB]: Public Health Bioinformatics is the GitHub repository that contains the TheiaProk workflows.

## Required Inputs

To activate `tbp-parser` you must set the following variable to true:

| Terra Task name | Variable | Type | Description | Default Value |
| :-------------- | :------- | :--- | :------------ | :---------- |
| `merlin_magic` | **call_tbp_parser** | Boolean | Set to `true` to activate `tbp-parser` | `false` |

## Optional Inputs

The following optional inputs are also available for user modification on Terra:

| Terra Task name | Variable | Type | Description | Default Value |
| :-------------- | :------- | :--- | :------------ | :---------- |
| `merlin_magic` | **tbp_parser_add_cs_lims** | Boolean | Set to `true` to add Cycloserine (CS) fields to the LIMS report | `false` |
| `merlin_magic` | **tbp_parser_coverage_regions_bed** | File | A BED file containing the regions to calculate percent coverage for | [tbdb-modified-regions.md](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb-modified-regions.bed) |
| `merlin_magic` | **tbp_parser_coverage_threshold** | Int | The minimum percentage of a region that has depth above the threshold set by `min_depth` (used for a gene/locus to pass QC) | 100 |
| `merlin_magic` | **tbp_parser_debug** | Boolean | Set to `false` to turn off debug mode for `tbp-parser` | `true` |
| `merlin_magic` | **tbp_parser_docker_image** | String | The Docker image to use when running `tbp-parser` | "us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:2.4.0" |
| `merlin_magic` | **tbp_parser_etha237_frequency** | Float | Minimum frequency for a mutation in ethA at protein position 237 to pass QC in `tbp-parser` | 0.1 |
| `merlin_magic` | **tbp_parser_expert_rule_regions_bed** | File | A file that contains the regions where R mutations and expert rules are applied |  |
| `merlin_magic` | **tbp_parser_min_depth** | Int | Minimum depth for a variant to pass QC in tbp_parser | 10 |
| `merlin_magic` | **tbp_parser_min_frequency** | Int | The minimum frequency for a mutation to pass QC | 0.1 |
| `merlin_magic` | **tbp_parser_min_read_support** | Int | The minimum read support for a mutation to pass QC | 10 |
| `merlin_magic` | **tbp_parser_operator** | String | Fills the "operator" field in the tbp_parser output files | "Operator not provided" |
| `merlin_magic` | **tbp_parser_output_seq_method_type** | String | Fills out the "seq_method" field in the tbp_parser output files | "Sequencing method not provided" |
| `merlin_magic` | **tbp_parser_rpob449_frequency** | Float | Minimum frequency for a mutation at protein position 449 to pass QC in `tbp-parser` | 0.1 |
| `merlin_magic` | **tbp_parser_rrl_frequency** | Float | Minimum frequency for a mutation in rrl to pass QC in `tbp-parser` | 0.1 |
| `merlin_magic` | **tbp_parser_rrl_read_support** | Int | Minimum read support for a mutation in rrl to pass QC in `tbp-parser` | 10 |
| `merlin_magic` | **tbp_parser_rrs_frequency** | Float | Minimum frequency for a mutation in rrs to pass QC in `tbp-parser` | 0.1 |
| `merlin_magic` | **tbp_parser_rrs_read_support** | Int | Minimum read support for a mutation in rrs to pass QC in `tbp-parser` | 10 |
| `merlin_magic` | **tbp_parser_tngs_data** | Boolean | Set to `true` to enable tNGS-specific parameters and runs in `tbp-parser` | `false` |

[Find the outputs for `tbp-parser` in TheiaProk on Terra here](../outputs/theiaprok.md).
