---
title: Looker Report
---

The Looker report is intended for use in Google's Looker Data Studio for dashboarding purposes. It offers a highly condensed version of the resistance calls (using the `looker_interpretation` field from the laboratorian report) for a quick summary of the sample’s drug resistance profile.

## Resistance columns

Each antimicrobial drug listed in the input `--coverage_bed` file will be represented by one column in the Looker report. The column name will be the antimicrobial drug name in lowercase. This can be customized by adding content to the `FIND_AND_REPLACE` dictionary in a configuration file, though please be aware that this will change every instance of that drug name in **all** output reports.

Drug interpretation severity is ranked as follows (from high to low): **R > R-Interim > U > S-Interim > S > WT > Insufficient Coverage > NA**.

| Column name | Explanation | Source |
| --- | --- | --- |
| antimicrobial_drug_name (e.g., rifampicin) | The highest `looker_interpretation` identified for that drug in the Laboratorian report | Determined by tbp-parser |

## Miscellaneous columns

These miscellaneous columns are also included in the Looker report to provide additional context about the sample. They can be customized by adding content to the `FIND_AND_REPLACE` dictionary in a configuration file.

| Column name | Explanation |
| --- | --- |
| sample_id | The name of the sample |
| output_seq_method_type | The sequencing method used to generate the data; can be set with the `--sequencing_method` input parameter. If left blank, "Sequencing method not provided" is the default value |
| lineage | The lineage of the sample (the `main_lin` field as reported by TBProfiler); for example, lineage1.2.1.2.1  |
| ID | The lineage of the sample in human-readable language (the same as the `Lineage ID` column in the LIMS report; see [the relevant documentation](./lims.md#lineage-id-language) for more information) |
| analysis_date | The date `tbp-parser` was run in YYYY-MM-DD HH:MM format |
| operator | The name of the person who ran `tbp-parser`; can be provided with the `--operator` input parameter. If left blank, "Operator not provided" is the default value. |

## Customizing column names

To overwrite any of the output column names or text in the laboratorian report, please use the following format in a [configuration file](../inputs.md#configuration-file) or use the command-line parameter `--find_and_replace`:

```yaml
FIND_AND_REPLACE:
  "sample_id": "My_Sample_ID_Column"
  "rifampicin": "rifampin"
```

Please note that this will rename every instance of that text in **all** output reports (every instance of "sample_id" will be renamed to "My_Sample_ID_Column" in all output files, etc.).
