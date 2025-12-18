---
title: LIMS Report
---

The LIMS report contains a summary of the mutations that confer resistance to any respective antimicrobial drug.

## Resistance columns

The LIMS report will contain one column for each drug defined in the `DRUG_COLUMNS_TO_GENE_COLUMNS` dictionary. Each drug will have one or more associated gene columns, depending on how many genes are responsible for resistance to that drug and indicated in the dictionary.

By default, **every** gene-drug combination found in the default `tbdb.bed` file is included in the LIMS report by virtue of its default inclusion in the `DRUG_COLUMNS_TO_GENE_COLUMNS` dictionary.

This dictionary takes the following format in `globals.py`:

```python
DRUG_COLUMNS_TO_GENE_COLUMNS = {
    "antimicrobial_drug_name": {
        "antimicrobial_column_name_in_lims_report": {
            "gene_name": "gene_drug_combo_column_name_in_lims_report",
            "gene2_name": "gene2_drug_combo_column_name_in_lims_report",            
            ...
        }
    },
    "antimicrobial2_drug_name": {
        "antimicrobial2_column_name_in_lims_report": {
            "gene3_name": "gene3_drug_combo_column_name_in_lims_report",
            ...
        }
    },
    ...
}
```

The column content is as follows:

| Column name  | Explanation | Source |
| --- | --- | --- |
| `antimicrobial_column_name_in_lims_report` | The highest `mdl_interpretation` identified for that drug in the Laboratorian report for the associated genes | Determined by tbp-parser |
| `gene_drug_combo_column_name_in_lims_report` | Any non-S (susceptible) mutations found in this gene that are responsible for the predicted resistance for the antimicrobial | Determined by tbp-parser |

These columns will be repeated for each drug and associated genes defined in the `DRUG_COLUMNS_TO_GENE_COLUMNS` dictionary. All genes associated with an antimicrobial drug will be grouped together with that antimicrobial column in the report.

### **Customizing the resistance column names**

In the example dictionary shown above, please do not modify these fields unless they are being removed or added as new drugs/genes:

- `antimicrobial_drug_name` is the name of the drug **as it appears in TBProfiler** (for example, "rifampicin").
- `gene_name` is the name of the gene **as it appears in TBProfiler** (for example, "rpoB").

However, the following fields can be customized to contain any text according to your laboratory's needs:

- `antimicrobial_column_name_in_lims_report` is the **desired name of the output column** in the LIMS report that indicates the highest resistance interpretation for that drug (for example, "RIF_RESISTANCE").
- `gene_drug_combo_column_name_in_lims_report` is the **desired name of the output column** in the LIMS report that indicates any mutations found in that gene that are responsible for the predicted resistance for that drug (for example, "RIF_RPOB_MUTATIONS").

To overwrite this dictionary in a configuration file, use the following format:

```yaml
DRUG_COLUMNS_TO_GENE_COLUMNS:
    rifampicin: # do not modify this text
        RIF_RESISTANCE: # this is the column name for the antimicrobial
            rpoB: RIF_RPOB_MUTATIONS # do not modify rpoB, but modify RIF_RPOB_MUTATIONS as needed
    linezolid:
        LZD_RESISTANCE:
            rrl: LZD_RRL_MUTATIONS
            rplC: LZD_RPLC_MUTATIONS
    ...
```

## Miscellaneous columns

These miscellaneous columns are not set with the dictionary described above and contain default names, but they can be renamed via the `OUTPUT_RENAMING` dictionary in a configuration file.

| Column name | Explanation | Source |
| --- | --- | --- |
| Sample_Name | The name of the sample | TBProfiler `"id"` field |
| Lineage_ID | The lineage of the sample in human-readable language | Determined by tbp-parser using the TBProfiler `"main_lineage"` and `"sub_lineage"` fields |
| Analysis_Date | The date `tbp-parser` was run in YYYY-MM-DD HH:SS format | Determined by tbp-parser at runtime |
| Operator | The name of the person who ran `tbp-parser`; can be provided with the `--operator` input parameter. If left blank, “Operator not provided” is the default value. | tbp-parser input parameter |
| Lineage | The lineage of the sample | TBProfiler `"main_lineage"` field |

### **Customizing miscellaneous column names**

To customize these column names in a configuration file, use the following format:

```yaml
OUTPUT_RENAMING:
  Sample_Name: "My_Sample_Column"
  Lineage_ID: "My_Lineage_ID_Column"
  Analysis_Date: "My_Analysis_Date_Column"
  Operator: "My_Operator_Column"
  Lineage: "My_Lineage_Column"
```

Please note that this will rename every instance of that text in all output reports (all "Sample_Name" will be renamed to "My_Sample_Column" in all output files, etc.).
