---
title: "Coverage Report"
---

The coverage report lists every gene and its percent gene coverage over a minimum depth (default: 10) relative to the H37Rv genome.

Please note that user-provided coverage regions always take precedence over default values.

## WGS Coverage Report

| Column name | Explanation |
| :---------- | :---------- |
| Gene | The name of the gene or locus |
| Percent_Coverage | The percent of the gene’s coding region that has a read depth over the minimum value (default: 10; user-customizable by altering `--min_depth`) |
| Warning | Indicates if any deletions were identified in the gene which may contribute to lower than expected coverage |

If run using the TheiaProk workflow series, there will be an additional column that contains only the name of the sample, which is useful when concatenating many reports as it helps differentiate which gene belongs to which sample.

## tNGS-specific information

If the `--tngs` flag is used, the report contains the following fields:

| Column name | Explanation |
| :---------- | :---------- |
| Gene | The name of the gene or locus |
| Coverage_Breadth_reportableQC_region | The percent of the gene (positions determined by the regions covered by the tNGS Deeplex + CDPH assay primers that are considered reportable by CDPH) that is covered at a depth greater than the `--min_depth` value |
| QC_Warning | Indicates if any deletions were identified in the gene which may contribute to lower than expected coverage |
| Coverage_Breadth_R_expert-rule_region | The percent of the regions (positions that could contain any resistance-conferring mutations or require expert-rule application) that is covered at a depth greater than the `--min_depth` value |

Coverage regions are determined with either the default [/data/tbdb-modified-regions.bed](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb-modified-regions.bed) (collected on Sep 1, 2023 from the TBProfiler repository, or if `--tngs`, [/data/tngs-reportable-regions.bed](https://github.com/theiagen/tbp-parser/blob/main/data/tngs-reportable-regions.bed).

The R-expert rule region is determined only if `--tngs` is indicated and uses the ranges in [/data/tbdb-expert-regions.bed](https://github.com/theiagen/tbp-parser/blob/main/data/tbdb-expert-regions.bed).
