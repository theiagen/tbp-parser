---
title: TheiaProk Inputs on Terra
---

When running `tbp-parser` as part of the TheiaProk workflow series ([find documentation for TheiaProk here](https://theiagen.notion.site/Theiagen-Public-Health-Resources-a4bd134b0c5c4fe39870e21029a30566?pvs=4)) on [Terra.bio](https://terra.bio), an optional input must be activated to instruct TheiaProk to run `tbp-parser`.

`tbp-parser` is not on by default due to the nature of this tool and its outputs.

!!! info annotate "TheiaProk Version"
    This information only corresponds to PHB v2.2.0. These inputs and outputs may not be applicable to other versions of TheiaProk.

*[PHB]: Public Health Bioinformatics is the GitHub repository that contains the TheiaProk workflows.

## Required Inputs

To activate `tbp-parser` you must set the following variable to true:

| Terra Task name | Variable | Type | Default value | Description |
| :-------------- | :------- | :--- | :------------ | :---------- |
| `merlin_magic` | `tbprofiler_additional_outputs` | Boolean | `false` | Set to `true` to activate `tbp-parser` |

## Optional Inputs

The following optional inputs are also available for user modification on Terra:

| Terra Task name | Variable | Type | Default value | Description |
| :-------------- | :------- | :--- | :------------ | :---------- |
| `merlin_magic` | `tbp_parser_output_seq_method_type` | String | "WGS" | Fills out the “seq_method” field in the tbp_parser output files |
| `merlin_magic` | `tbp_parser_operator` | String | "Operator not provided" | The operator who ran the analysis; used in the LIMS & Looker reports |
| `merlin_magic` | `tbp_parser_min_depth` | Int | 10 | The minimum depth of coverage required for a site to pass QC |
| `merlin_magic` | `tbp_parser_min_frequency` | Int | 0.1 | The minimum frequency for a mutation to pass QC (0.1 -> 10%) |
| `merlin_magic` | `tbp_parser_min_read_support` | Int | 10 | The minimum read support for a mutation to pass QC |
| `merlin_magic` | `tbp_parser_coverage_threshold` | Int | 100 | The minimum percentage of a region that has depth above the threshold set by `min_depth` (used for a gene/locus to pass QC) |
| `merlin_magic` | `tbp_parser_coverage_regions_bed` | File | [tbdb-modified-regions.md](https://github.com/theiagen/tbp-parser/blob/v1.6.0/data/tbdb-modified-regions.bed) | A BED file containing the regions to calculate percent coverage for |
| `merlin_magic` | `tbp_parser_debug` | Boolean | false | Turn on debug mode for tbp-parser |
| `merlin_magic` | `tbp_parser_add_cs_lims` | Boolean | false | Adds Cycloserine (CS) fields to the LIMS report |
| `merlin_magic` | `tbp_parser_docker_image` | String | "us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:1.6.0" | The Docker image to use when running tbp-parser |

[Find the outputs for `tbp-parser` in TheiaProk on Terra here](../outputs/theiaprok.md).