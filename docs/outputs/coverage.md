---
title: "Coverage Report"
---

The coverage report lists every gene and its percent gene coverage over a minimum depth (default: 10) relative to the H37Rv genome; regions are determined by using the `--coverage_bed` input file. This report is useful for determining whether any genes of interest have low coverage that may impact the reliability of resistance calls.

## Coverage report columns

| Column name | Explanation |
| :---------- | :---------- |
| sample_name | The name of the sample |
| locus_tag | The locus tag for the region |
| gene_name | The name of the gene |
| percent_coverage | The percentage of a region (specified by the `--coverage_bed` input file) that has a read depth over the minimum value (default: 10; user-customizable by altering `--min_depth`) |
| average_depth | The average read depth across that region (specified by the `--coverage_bed` input file) |
| err_percent_coverage | The percentage of the region (specified by the `--err_coverage_bed` input file) that has a read depth over the minimum value (default: 10; user-customizable by altering `--min_depth`)<br>_Only appears if `--err_coverage_bed` is provided_ |
| err_average_depth | The average read depth across a region (specified by the `--err_coverage_bed` input file)<br>_Only appears if `--err_coverage_bed` is provided_ |
| qc_warning | Indicates if any deletions were identified in the gene which may contribute to lower than expected coverage |

### Locus vs Target coverage reports

Two coverage reports are generated: 

- a **locus** coverage report, which aggregates multiple target regions if they belong to the same locus tag if neccessary (for example, during a tNGS analysis, several genes are covered by mulitple primers)
- a **target** coverage report, which keeps each record in the provided `--coverage_bed` file separate (for example, during a tNGS analysis, each primer region is kept separate)

If each locus is covered by only one target region, the two reports will be identical. 


### Example Coverage Report

In the following example, katG is covered by four different target regions (katG1, katG2, katG3, and katG4) in the `--coverage_bed` file that overlap. In the locus coverage report, these four regions are aggregated into one record for the katG gene, while in the target coverage report, each of the four regions is kept separate.

For regions that do not overlap, such as rrl and rpoB, the locus report only looks at the individual regions sequenced and not the area in between.

Here is an excerpt from a locus coveage report from a tNGS analysis:

```text
sample_name	locus_tag	gene_name	percent_coverage	average_depth	err_percent_coverage	err_average_depth	qc_warning
sample01	Rv1305	atpE	100	1882.364	100	2409.771	
sample01	Rv2416c	eis	100	2290.005	100	2538.581	
sample01	Rv3795	embB	100	2376.114	100	2919.096	
sample01	Rv0006	gyrA	100	2313.786	100	2823.435	
sample01	Rv0005	gyrB	100	1583.527	100	1949.587	
sample01	Rv1484	inhA	100	2326.662	100	2551.08	
sample01	Rv1908c	katG	100.000	2194.408	100.000	2333.477	
sample01	Rv0678	mmpR5	100	3158.585	100	3674.649	
sample01	Rv2043c	pncA	100	1715.247	100	1829.639	
sample01	Rv0701	rplC	100	2889.414	100	3095.667	
sample01	Rv0667	rpoB	100	3285.271	100	2342.436	Contains valid deletion(s): c.1291_1292delAGinsCC; c.1291_1292delAGinsGC; c.1291_1292delAGinsTT; c.1327_1329delTTGinsCCC; c.1291_1292delAGinsTC; c.1327_1329delTTGinsCTC; c.1327_1329delTTGinsCTT
sample01	EBG00000313339	rrl	100	2383.529	100	2952.955
sample01	EBG00000313325	rrs	100	3128.388	100	3719.095	
```

Here is a corresponding exerpt from a target coverage report:

```text
sample_name	locus_tag	gene_name	percent_coverage	average_depth	err_percent_coverage	err_average_depth	qc_warning
sample01	Rv1305	atpE	100	1882.364	100	2409.771	
sample01	Rv2416c	eis	100	2290.005	100	2538.581	
sample01	Rv3795	embB	100	2376.114	100	2919.096	
sample01	Rv0006	gyrA	100	2313.786	100	2823.435	
sample01	Rv0005	gyrB	100	1583.527	100	1949.587	
sample01	Rv1484	inhA	100	2326.662	100	2551.08	
sample01		katG1	100	1840.63	100	2032.022	
sample01		katG2	100	1941.742	100	2579.653	
sample01		katG3	100	1913.291	100	2449.762	
sample01		katG4	100	1922.261	100	2353.081	
sample01	Rv0678	mmpR5	100	3158.585	100	3674.649	
sample01	Rv2043c	pncA	100	1715.247	100	1829.639	
sample01	Rv0701	rplC	100	2889.414	100	3095.667	
sample01		rpoB_1	100	4851.581	100	6403.333
sample01		rpoB_2	100	1901.996	100	2281.826
sample01		rrl_1	100	1865.669	100	2085.516	
sample01		rrl_2	100	2696.711	100	3166.373	
sample01	EBG00000313325	rrs	100	3128.388	100	3719.095	
```

## Customizing column names

To overwrite any of the column names in a configuration file, use the following format:

```yaml
FIND_AND_REPLACE:
  "gene_name": "My_Gene_Column"
  "percent_coverage": "My_Percent_Coverage_Column"
  ...
```

Please note that this will rename every instance of that text in **all** output reports (every instance of "gene_name" will be renamed to "My_Gene_Column" in all output files, etc.).
