---
title: "Coverage Report"
---

The coverage report lists every gene and its percent gene coverage over a minimum depth (default: 10) relative to the H37Rv genome; regions are determined by using the `--tbdb_bed` input file. This report is useful for determining whether any genes of interest have low coverage that may impact the reliability of resistance calls.

## WGS Coverage Report

| Column name | Explanation |
| :---------- | :---------- |
| Gene | The name of the gene or locus |
| Percent_Coverage | The percent of the gene’s coding region that has a read depth over the minimum value (default: 10; user-customizable by altering `--min_depth`) |
| Average_Locus_Coverage | The average read depth across the entire gene |
| Warning | Indicates if any deletions were identified in the gene which may contribute to lower than expected coverage |

## tNGS-specific information

If the `--tngs` flag is used, the report contains the following fields:

| Column name | Explanation |
| :---------- | :---------- |
| Gene | The name of the gene or locus |
| Coverage_Breadth_reportableQC_region | The percent of the gene that is covered at a depth greater than the `--min_depth` value; only for the region specified in the `--tbdb_bed` |
| Average_Locus_Coverage | The average read depth across the entire gene |
| QC_Warning | Indicates if any deletions were identified in the gene which may contribute to lower than expected coverage |
| Coverage_Breadth_R_expert-rule_region | The percent of the gene that is covered at a depth greater than the `--min_depth` value but only for the region of the gene indicated by the `--err_bed` file |

## Customizing column names

To overwrite any of the column names in a configuration file, use the following format:

```yaml
FIND_AND_REPLACE:
  "Gene": "My_Gene_Column"
  "Percent_Coverage": "My_Percent_Coverage_Column"
  ...
```

Please note that this will rename every instance of that text in all output reports (every instance of "Gene" will be renamed to "My_Gene_Column" in all output files, etc.).
