---
title: Looker Report
---

The Looker report is intended for use in Google's Looker Data Studio for dashboarding purposes. It offers a highly condensed version of the resistance calls (using the `looker_interpretation` field from the laboratorian report) for a quick summary of the sample’s drug resistance profile.

## Resistance columns

Each antimicrobial drug listed in the input `--tbdb_bed` file will be represented by one column in the Looker report. The column name will be the antimicrobial drug name in lowercase. This can be customized by adding content to the `OUTPUT_RENAMING` dictionary in a configuration file, though please be aware that this will change every instance of that drug name in all output reports.

| Column name | Explanation | Source |
| --- | --- | --- |
| antimicrobial_drug_name (e.g., rifampicin) | The highest `looker_interpretation` identified for that drug in the Laboratorian report | Determined by tbp-parser |

## Miscellaneous columns

These miscellaneous columns are also included in the Looker report to provide additional context about the sample. They can be customized by adding content to the `OUTPUT_RENAMING` dictionary in a configuration file.

| Column name | Explanation |
| --- | --- |
| sample_id | The name of the sample |
| output_seq_method_type | The sequencing method used to generate the data; can be set with the `--sequencing_method` input parameter. If left blank, “Sequencing method not provided” is the default value |
| lineage | The lineage of the sample (the `main_lin` field as reported by TBProfiler); for example, lineage1.2.1.2.1  |
| ID | The lineage of the sample in human-readable language (the same as `M_DST_A01_ID` in the LIMS report) |
| analysis_date | The date `tbp-parser` was run in YYYY-MM-DD HH:SS format |
| operator | The name of the person who ran `tbp-parser`; can be provided with the `--operator` input parameter. If left blank, “Operator not provided” is the default value. |

## Customizing column names

To customize any of the column names, you must use the following format in a configuration file:

```yaml
OUTPUT_RENAMING:
  sample_id: "My_Sample_ID_Column"
  rifampicin: "rifampin"
```

This will rename the `sample_id` column to `My_Sample_ID_Column` and the `rifampicin` column to `rifampin` in the Looker report, but also rename every instance of that text in all output reports (all "rifampicin" will be renamed to "rifampin").
