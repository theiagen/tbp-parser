---
title: LIMS Report
---

The LIMS report contains a summary of the mutations that confer resistance to any respective antimicrobial drug.

The report can be divided into two sections: resistance columns and miscellaneous columns. Resistance columns are determined by the `lims_report_format_yml` file, which defines the antimicrobial drugs and associated genes that should be included in the report. Miscellaneous columns include more generic information such as the sample name, lineage, analysis date, and operator.

## Resistance columns

The LIMS report will contain one column for each drug defined by the `lims_report_format_yml` file. Each drug will have one or more associated gene columns, depending on how many genes are responsible for resistance to that drug.

Drug interpretation severity is ranked as follows (from high to low): **R > Insufficient Coverage > U > S > WT**. 

| Column names  | Explanation | Source |
| --- | --- | --- |
| `antimicrobial_column_name_in_lims_report` | The highest `mdl_interpretation` identified for that drug in the Laboratorian report for the associated genes | Determined by tbp-parser |
| `gene_drug_combo_column_name_in_lims_report` | Any non-S (susceptible) mutations found in this gene that are responsible for the predicted resistance for the antimicrobial | Determined by tbp-parser |

These columns will be repeated for each drug and associated genes defined in the `lims_report_format.yml` file. All genes associated with an antimicrobial drug will be grouped together with that antimicrobial's column in the report.

By default, **every** gene-drug combination found in the default `tbdb.bed` file is included in the LIMS report by virtue of its default inclusion in the `lims_report_format_yml` file.

### **Syntax and logic used**

The language used for the resistance columns is the same between tNGS and WGS analysis. All mutations that failed quality in the position are not considered for the LIMS report.

Please note that drug interpretation severity is ranked as follows (from high to low): **R > Insufficient Coverage > U > S > WT**. 

| Column Type | Language | Explanation | Responsible MDL Interpretation(s) |
|---|---|---|---|
| Drug | Mutations(s) associated with resistance to <antimicrobial> detected | The highest severity mutation for any of the gene(s) associated with that antimicrobial drug has an "R" (resistant) MDL interpretation in the Laboratorian report | R |
| Drug | The detected mutation(s) have uncertain significance. Resistance to <antimicrobial> cannot be ruled out | The highest severity mutation found any of the gene(s) associated with that antimicrobial drug has a "U" (uncertain significance) MDL interpretation in the Laboratorian report | U |
| Gene-drug combo | p.Ala689Val (e.g) | A list of the mutations for that gene-drug combination with either an "R" or "U" MDL interpretation in the Laboratorian report. The format of the mutation is "p.(amino acid change)" by default or "n.(nucleotide change)" if no protein change occurred. | R, U |
| Drug | No mutations associated with resistance to <antimicrobial> detected | The highest severity mutation for any of the gene(s) associated with that antimicrobial drug has an "S" (susceptible) or "WT" (wild-type) MDL interpretation in the Laboratorian report | S, WT |
| Gene-drug combo | No high confidence mutations detected | The highest severity mutation for that gene-drug combination has an "S" MDL interpretation. | S |
| Gene-drug combo | No mutations detected | The highest severity mutation for that gene-drug combination has a "WT" or "NA" MDL interpretation; no mutations were found in that gene, or any mutations failed positional QC. | WT |
| Drug | Pending Retest | The highest severity mutation for any of the gene(s) associated with that antimicrobial drug has an "Insufficient Coverage" MDL interpretation in the Laboratorian report | Insufficient Coverage |
| Gene-drug combo | No sequence | The highest severity mutation for that gene-drug combination has an "Insufficient Coverage" MDL interpretation | Insufficient Coverage |

In tbp-parser, rifampicin uses slightly different language:

| Column Type | Language | Explanation | Responsible MDL Interpretation |
|---|---|---|---|
| Drug | Predicted low-level resistance to rifampicin. May test susceptible by phenotypic methods | One or more of the following mutations in rpoB were detected: Leu430Pro, Asp435Tyr, His445Asn, His445Ser, His445Leu, His445Cys, Leu452Pro, Ile491Phe. No other R mutations were found in rpoB. | R |
| Drug | Predicted resistance to rifampicin | An R mutation was detected in rpoB that is not part of the list above. | R |
| Drug | Predicted susceptibility to rifampicin. The detected synonymous mutation(s) do not confer resistance | The highest severity mutation for the gene(s) associated with rifampicin has an "S" MDL interpretation and synonymous mutations were present within the rpoB RRDR region. | S |
| Drug | Predicted susceptibility to rifampicin | The highest severity mutation for the gene(s) associated with rifampicin has an "S" MDL interpretation and **no** synonymous mutations were detected in the rpoB RRDR region. | S |
| Gene-drug combo | p.Ala689Ala [synonymous] (e.g.) | If a [synonymous] tag follows a mutation, it is from the rpoB RRDR region. | S |

### **Customizing resistance column names**

The output column names can be customized to contain any text according to your laboratory's needs by providing a custom `lims_report_format_yml` file, which should take the following format:

```yaml
# do not modify unbracketed text
# <this text can be fully customized>
# [this text must match TBProfiler nomenclature for drug and gene names]

- drug: [drug_name]
  drug_code: <antimicrobial_column_name_in_lims_report>
  gene_codes:
    [gene_name]: <column_name_for_gene_drug_combo_in_lims_report> 
    [gene_name]: <column_name_for_gene_drug_combo_in_lims_report>
    ...
- drug: [drug_name]
  drug_code: <antimicrobial_column_name_in_lims_report>
  gene_codes: {}
...
```

- `drug_name` is the name of the drug **as it appears in TBProfiler** (for example, "rifampicin").
- `gene_name` is the name of the gene **as it appears in TBProfiler** (for example, "rpoB").
- `antimicrobial_column_name_in_lims_report` is the **desired name of the output column** in the LIMS report that indicates the highest resistance interpretation for that drug (for example, "RIF").
- `column_name_for_gene_drug_combo_in_lims_report` is the **desired name of the output column** in the LIMS report that indicates any mutations found in that gene that are responsible for the predicted resistance for that drug (for example, "RIF_rpoB").

For example:

```yaml
- drug: rifampicin
  drug_code: RIF
  gene_codes:
    rpoB: RIF_rpoB
- drug: amikacin
  drug_code: AMK
  gene_codes:
    bacA: AMK_bacA
    ccsA: AMK_ccsA
    eis: AMK_eis
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

### **Lineage ID language**

Lineage ID is only designated if the percentage of LIMS genes that pass locus QC is greater than the minimum percentage, indicated by `--min_percent_loci_covered` (default 0.7 -> 70%).

The language used is different between tNGS and WGS. 

| WGS Language | Explanation |
| --- | --- | 
| DNA of Mycobacterium tuberculosis species detected | The TBProfiler `"main_lineage"` field contains "lineage" |
| DNA of Mycobacterium bovis BCG detected | The TBProfiler `"main_lineage"` field or `"sub_lineage"` field contains "BCG" |
| DNA of Mycobacterium bovis (not BCG) detected | The TBProfiler `"main_lineage"` field or `"sub_lineage"` field does **not** contain "BCG" but _does_ contain "bovis" or "La1" |
| DNA of Mycobacterium tuberculosis complex detected | The TBProfiler `"main_lineage"` field is blank, "NA", or does not exist |
| DNA of Mycobacterium tuberculosis NOT detected | The `--min_percent_loci_covered` QC check fails |

| tNGS Language | Explanation |
|---|---|
| DNA of Mycobacterium bovis BCG detected | pncA contains a "p.His57Asp" mutation that passes positional QC |
| DNA of Mycobacterium tuberculosis complex detected (M. bovis not ruled out) | pncA contains a "p.His57Asp" mutation that **does not** pass positional QC |
| DNA of Mycobacterium tuberculosis complex detected (not M. bovis) | pncA **does not** have a "p.His57Asp" mutation |
| DNA of Mycobacterium tuberculosis NOT detected | The `--min_percent_loci_covered` QC check fails |

### **Customizing miscellaneous column names**

To overwrite any of the output column names or text in the laboratorian report, please use the following format in a [configuration file](../inputs.md#configuration-file) or use the command-line parameter `--find_and_replace`:

```yaml
FIND_AND_REPLACE:
  rifampicin: "rifampin"
  fbiD: "Rv2983"
  mmpR5: "Rv0678"
  "Sample Name": "sample accession"
  "DNA of Mycobacterium bovis BCG detected" : "M. bovis BCG"
```

Please note that this will rename every instance of that text in **all** output reports (every instance of "Sample Name" will be renamed to "sample accession" in all output files, etc.).

## Example LIMS report

The following is an example of the default LIMS report (transposed for readability).

```text
Sample Name,sample01
Sample ID,DNA of Mycobacterium tuberculosis species detected
AMK,Mutation(s) associated with resistance to amikacin detected
AMK_bacA,No sequence
AMK_ccsA,No high confidence mutations detected
AMK_eis,No mutations detected
AMK_rrs,n.1401A>G
PZA,The detected mutation(s) have uncertain significance. Resistance to pyrazinamide cannot be ruled out
PZA_clpC1,No mutations detected
PZA_rpsA,No high confidence mutations detected
PZA_Rv1258c,p.Glu194fs
RIF,Pending Retest
RIF_rpoB,No high confidence mutations detected
RIF_rpoC,No mutations detected
RIF_Rv1129c,No sequence
Analysis date,2026-03-05 19:36
Operator,Operator not provided
Lineage,lineage2
```
