---
title: Looker Report
---

The Looker report is intended for use in Google's Looker Studio Data Studio for dashboarding purposes. It offers a highly condensed version of the resistance calls (using the `looker_interpretation` field from the laboratorian report) for a quick summary of the sample’s drug resistance profile.

### Explanation of column headers

| Column name | Explanation |
| --- | --- |
| sample_id | The name of the sample |
| output_seq_method_type | The sequencing method used to generate the data; can be set with the `--sequencing_method` input parameter. If left blank, “Sequencing method not provided” is the default value |
| amikacin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| bedaquiline | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| capreomycin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| clofazimine | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| ethambutol | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| ethionamide | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| isoniazid | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| kanamycin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| levofloxacin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| linezolid | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| moxifloxacin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| pyrazinamide | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| rifampin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| streptomycin | The highest `looker_interpretation` resistance identified for mutations associated with this drug |
| lineage | The lineage of the sample (the `main_lin` field as reported by TB-Profiler); for example, lineage1.2.1.2.1  |
| ID | The lineage of the sample in human-readable language (the same as `M_DST_A01_ID` in the LIMS report) |
| analysis_date | The date `tbp-parser` was run in YYYY-MM-DD HH:SS format |
| operator | The name of the person who ran `tbp-parser`; can be provided with the `--operator` input parameter. If left blank, “Operator not provided” is the default value. |

Please note that occasionally, the `looker_interpretation` field can differ from the `mdl_interpretation` field. Typically, they are identical, but occasionally, the `mdl_interpretation` column will call a variant-drug combination “susceptible” (S), while the `looker_interpretation` column will call the same combination “uncertain” (U). Be aware of this difference when choosing an interpretation to report.

*[highest]: This is either resistant (R), uncertain (U), or susceptible (S)
