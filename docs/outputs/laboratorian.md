---
title: Laboratorian Report
---

The laboratorian report is the main report produced by `tbp-parser` and is used to generate all of the other reports. What follows is an explanation of all the columns in the report.

Any fields from TBProfiler are from the input JSON file produced by TBProfiler.

| Column name | Explanation | Column source |
| --- | --- | --- |
| sample_id | The name of the sample | TBProfiler "id" field |
| tbprofiler_gene_name | The name of the gene where the mutation has been identified | TBProfiler "gene_name" field |
| tbprofiler_locus_tag | The locus tag for the mutation that has been identified | TBProfiler "locus_tag" field _OR_ the GENE_TO_LOCUS_TAG dictionary created from the --tbdb_bed input file if the field is blank |
| tbprofiler_variant_substitution_type | The type of mutation identified, whether or not it was a frameshift, missense, or synonymous mutation | TBProfiler "type" field |
| tbprofiler_variant_substitution_nt | The mutation in nucleotide format | TBProfiler "nucleotide_change" field |
| tbprofiler_variant_substitution_aa | The mutation in amino acid format, if possible | TBProfiler "protein_change" field |
| confidence | Contains either:<br>- the WHO annotation<br>- an indication that there is no WHO annotation<br>- NA for when there is no mutation | Edited by tbp-parser, originates from the TBProfiler "confidence" field |
| antimicrobial | The antimicrobial drug that may be affected by this mutation |  TBProfiler "annotation.drug" field, split into multiple rows if multiple annotation items are present. May also originate from the "gene_associated_drugs" field if not all are included in the annotation |
| looker_interpretation | The drug resistance interpretation intended for the Looker report | Determined by tbp-parser |
| mdl_interpretation | The drug resistance interpretation intended for the LIMS report | Determined by tbp-parser |
| depth | The depth of coverage at the mutation | TBProfiler "depth" field |
| frequency | The frequency of the mutation in the reads | TBProfiler "freq" field |
| read_support | How many reads support the mutation (depth * frequency) | Calculated by tbp-parser |
| rationale | Contains an indication of what was used (the WHO annotation, the specific expert rule used, or neither) to create the two interpretations | Determined by tbp-parser |
| warning | Any potential quality warnings that may indicate lower reliability | Determined by tbp-parser |
| gene_tier | The gene tier of the mutation’s gene (Tier 1, Tier 2, or NA) | Determined by the GENE_TO_TIER dictionary created from the --gene_tier_tsv input file |
| source | The source of the mutation information (WHO v2 catalogue, tbdb, etc.) | TBProfiler "annotation.source" field |
| tbdb_comment | Any comments from TBProfiler about the mutation | TBProfiler "annotation.comment" field |

Because of how a particular mutation may contribute resistance to different drugs at the same time, each mutation is listed multiple times, once for each antimicrobial drug that could be affected.

Any genes that do _not_ have any mutations are also included in the laboratorian report with NA or WT in the appropriate field. This results in a report with many rows and often, rows with very similar values.

The laboratorian report contains the "complete picture" of the sample and is incredibly useful for understanding the sample’s drug resistance profile.

## Customizing column names

To overwrite any of the column names in a configuration file, use the following format:

```yaml
OUTPUT_RENAMING:
  sample_id: "My_Sample_ID_Column"
  tbprofiler_gene_name: "My_Gene_Name_Column"
  ...
```

Please note that this will rename every instance of that text in all output reports (all "sample_id" will be renamed to "My_Sample_ID_Column" in all output files, etc.).
