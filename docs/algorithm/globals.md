---
title: Global Variables
---

The following global variables are used in the `tbp-parser` algorithm. These variables are used to determine columns in the output files, thresholds for quality control, and various other look-up tables for the algorithm.

These variables can be overwritten using a configuration file. For more information on how to overwrite these variables, please see the [Configuration File](../inputs/inputs.md#configuration-file) section.

## Global Variables

The content of the global variables can be found by browsing the [globals.py](https://github.com/theiagen/tbp-parser/blob/main/tbp_parser/globals.py) file on GitHub.

Some global variables are excluded as their modifications will either have no effect or will break the algorithm without further changes.

| Variable Name | Description | Format |
| :------------ | :---------- | :------------ |
| ANNOTATION_TO_INTERPRETATION | A dictionary to turn TBProfiler WHO annotations into their corresponding CDPH Looker or MDL interpretations | Dict[Dict[String: String]] |
| ANTIMICROBIAL_CODE_TO_DRUG_NAME | A dictionary that matches CDPH LIMS antimicrobial codes to the corresponding drug name; used to create the LIMS report | Dict[String: String] |
| ANTIMICROBIAL_CODE_TO_DRUG_NAME_CS | A dictionary that matches the LIMS antimicrobial code to the relevant antimicrobial drug name for cycloserine; activated by the `--add_cs_lims` flag | Dict[String: String] |
| DRUG_COLUMNS_TO_GENE_COLUMNS | The dictionary that matches CDPH LIMS antimicrobial codes to the corresponding genes and their corresponding CDPH LIMS codes; used to create the LIMS report | Dict[] |
| DRUG_COLUMNS_TO_GENE_COLUMNS_tNGS | A dictionary that matches CDPH LIMS antimicrobial codes to the corresponding genes and their corresponding CDPH LIMS codes for tNGS data; used to create the LIMS report; fills the `DRUG_COLUMNS_TO_GENE_COLUMNS` | Dict[Dict[String: String]] |
| DRUG_COLUMNS_TO_GENE_COLUMNS_WGS | A dictionary that matches CDPH LIMS antimicrobial codes to the corresponding genes and their corresponding CDPH LIMS codes for WGS data; used to create the LIMS report; fills the `DRUG_COLUMNS_TO_GENE_COLUMNS` | Dict[Dict[String: String]] |
| DRUG_COLUMNS_TO_GENE_COLUMNS_CS | A dictionary that adds CDPH LIMS antimicrobial codes to the corresponding genes and their corresponding CDPH LIMS codes for cycloserine; activated by the --add_cs_lims flag; used to create the LIMS report; is added to `DRUG_COLUMNS_TO_GENE_COLUMNS` | Dict[Dict[String: String]] |
| ANTIMICROBIAL_DRUG_NAME_LIST | A list of the antimicrobial drug names; used to create the Looker report | List[String] |
| ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME | A dictionary that matches the antimicrobial drugs to the genes that may confer resistance to them; see also <https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv>; used to create the Looker report | Dict[String: List[String]] |
| COVERAGE_THRESHOLD | The coverage threshold (the minimum percent coverage of a gene/locus over the minimum depth); See also MIN_DEPTH | Float |
| ETHA237_FREQUENCY | The minimum frequency for a mutation to be considered for ethA at protein position 237 | Float |
| GENE_LIST | A list of genes that correspond to a certain set of expert rules; Rv0678 is equivalent to mmpR5 | List[String] |
| GENE_TO_ANTIMICROBIAL_DRUG_NAME | A dictionary corresponding each gene to the drug they may confer resistance to, including the genes in the TBDB watchlist. See also: <https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv>, and <https://github.com/jodyphelan/tbdb/blob/master/tbdb.watchlist.csv>; used to create the Laboratorian report | Dict[String: List[String]] |
| GENE_TO_LOCUS_TAG | A dictionary that matches each gene to its corresponding locus tag; see also: <https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb.bed>; used to create the Laboratorian report | Dict[String: String] |
| GENE_TO_TIER | A dictionary that matches each gene to it's corresponding tier; tier information provided by CDPH; used to create the Laboratorian report | Dict[String: String] |
| GENES_FOR_LIMS | The list of genes used to generate the CDPH LIMS report; used to create the LIMS report | List[String] |
| GENES_FOR_LIMS_tNGS | The list of genes used to generate the CDPH LIMS report for tNGS data; used to create the LIMS report; fills the `GENES_FOR_LIMS` | List[String] |
| GENES_FOR_LIMS_WGS | The list of genes used to generate the CDPH LIMS report for WGS data; used to create the LIMS report; fills the `GENES_FOR_LIMS` | List[String] |
| GENES_FOR_LIMS_CS | The list of genes used to generate the CDPH LIMS report for cycloserine; activated by the --add_cs_lims flag; used to create the LIMS report; is added to `GENES_FOR_LIMS` | List[String] |
| MIN_DEPTH | The minimum depth of a gene/locus to pass QC | Int |
| MIN_FREQUENCY | The minimum frequency of a mutation to pass QC | Float |
| MIN_READ_SUPPORT | The minimum read support of a mutation to pass QC | Int |
| OPERATOR | Used to fill the `OPERATOR` column in the output files | String |
| WHOV2_PROMOTER_REGIONS | A dictionary that contains regions that require special consideration | Dict[String: List[Int] _or_ List[List[Int]]] |
| RPOB449_FREQUENCY | The minimum frequency for a mutation to be considered for rpoB at protein position 449 | Float |
| RRL_FREQUENCY | The minimum frequency for a mutation to be considered for rrl | Float |
| RRL_READ_SUPPORT | The minimum read support for a mutation to be considered for rrl | Int |
| RRS_FREQUENCY | The minimum frequency for a mutation to be considered for rrs | Float |
| RRS_READ_SUPPORT | The minimum read support for a mutation to be considered for rrs | Int |
| SAMPLE_NAME | Used to fill the sample name column in the output files | String |
| SEQUENCING_METHOD | Used to fill the sequencing_method column in the output files | String |
| TNGS_REGIONS | The speciifc primer regions for the gene(s) that is(are) split into multiple sections | Dict[Dict[String: List[Int]]] |
