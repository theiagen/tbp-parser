---
title: Laboratorian Report
---

The laboratorian report is the main report produced by `tbp-parser` and is used to generate all of the other reports. What follows is an explanation of all the columns in the report.

### Explanation of column headers

| Column name | Explanation |
| --- | --- |
| sample_id | The name of the sample |
| tbprofiler_gene_name | The name of the gene where the mutation has been identified |
| tbprofiler_locus_tag | The locus tag for the mutation that has been identified |
| tbprofiler_variant_substitution_type | The type of mutation identified, whether or not it was a frameshift, missense, or synonymous mutation |
| tbprofiler_variant_substitution_nt | The mutation in nucleotide format |
| tbprofiler_variant_substitution_aa | The mutation in amino acid format, if possible |
| confidence | Contains either:<br>- the WHO annotation<br>- an indication that there was no WHO annotation<br>- NA for when there is no mutation |
| antimicrobial | The antimicrobial drug that may be affected by this mutation |
| looker_interpretation | The drug resistance interpretation intended for the Looker report |
| mdl_interpretation | The drug resistance interpretation intended for the LIMS report |
| depth | The depth of coverage at the mutation |
| frequency | The frequency of the mutation in the reads |
| read_support | How many reads support the mutation (depth * frequency) |
| rationale | Contains an indication of what was used (the WHO annotation, the specific expert rule used, or neither) to create the two interpretations |
| warning | Any potential quality warnings that may indicate lower reliability |
| gene_tier | The gene tier of the mutation’s gene (Tier 1, Tier 2, or NA) |

Because of how a particular mutation may contribute resistance to different drugs at the same time, each mutation is listed multiple times, once for each antimicrobial drug that could be affected. In addition, any genes that do *not* have any mutations are also included in the laboratorian report with NA or WT in the appropriate field. This results in a report with many rows and often, rows with very similar values. However, the laboratorian report contains the “complete picture” of the sample and is incredibly useful for understanding the sample’s drug resistance profile.
