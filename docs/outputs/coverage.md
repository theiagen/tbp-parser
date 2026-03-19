---
title: "Coverage Report"
---

The coverage report lists every gene and its percent gene coverage over a minimum depth (default: 10) relative to the H37Rv genome; regions are determined by using the `--coverage_bed` input file. This report is useful for determining whether any genes of interest have low coverage that may impact the reliability of resistance calls.

## Coverage report columns

ERR columns only appear if an `--err_coverage_bed` input file is provided. Please see the [Inputs](../inputs.md) for more details.

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

## Locus vs Target coverage reports

Each line in your BED file defines a **target region** — a stretch of the genome where tbp-parser counts reads and calculates coverage. How you define these regions directly controls what gets measured and what appears in your coverage reports.

- The **locus** coverage report aggregates all target regions that share the same locus tag and considers them as a single unit.
- The **target** coverage report keeps each BED file line as its own separate unit.

!!! info "Which report gets generated?"
    - If all BED file lines cover only one target region per locus, the two reports would be identical. In this case, only the locus coverage report is generated.
    - If any BED file lines cover multiple target regions per locus, both reports will be generated.

!!! info "QC applies locus-level coverage"
    All coverage-based QC decisions use the **locus-level** breadth of coverage. If katG1 has poor coverage but katG2 has sufficient coverage, the aggregated coverage may still pass or fail specific thresholds depending on the overall ratio.

### Simple case: one target region in BED file

```
BED file:
chrom    start    end    locus_tag    gene_name
Chrom    1        500    Rv0005       gyrB
```

Only the locus coverage report is produced:

```
locus_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv0005     gyrB       100.000           542.310

  → percent_coverage = (# positions ≥ min_depth) / 500
  → average_depth = (sum of reads at every position) / 500
```

### Split regions

A gene can be covered by multiple non-overlapping target regions (common in tNGS), resulting in multiple BED lines with the **same locus tag** but **different gene names** and a gap between them:

```
BED file:
chrom    start    end    locus_tag    gene_name
Chrom    100      200    Rv0667       rpoB_1
Chrom    300      400    Rv0667       rpoB_2

├── rpoB_1 ──┤              ├── rpoB_2 ──┤
100          200            300          400
               (gap: 201–299)
```

Because there are more targets (2) than loci (1), **both** reports are generated:

```
locus_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv0667     rpoB       95.000            280.500
                                   ↑ aggregated across positions 100–200 AND 300–400 (202 total positions)
```
```
target_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv0667     rpoB_1     100.000           350.100
sample01     Rv0667     rpoB_2     90.000            210.900
                                   ↑ each BED line measured independently
```

The gap (positions 201–299) is **not measured at all**, so no reads are counted there. The locus report aggregates only the positions defined by the BED lines (100–200 and 300–400), not the gap between them. This means the locus-level breadth of coverage reflects only the regions defined in the BED file, not the entire gene.

### Overlapping regions

A gene can be covered by multiple overlapping target regions (common in tNGS), resulting in multiple BED lines with the **same locus tag** but **different gene names** that overlap:

```
BED file:
chrom    start    end    locus_tag    gene_name
Chrom    1        200    Rv1908c      katG1
Chrom    150      400    Rv1908c      katG2

├──── katG1 ────┤
1               200
            ├────── katG2 ──────┤
            150                 400
```

Because there are now more targets (2) than loci (1), **both** reports are generated:

```
locus_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv1908c    katG       100.000            3471.463
                                   ↑ aggregated across ALL positions 1–400
```
```
target_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv1908c    katG1      100.000           2202.284
sample01     Rv1908c    katG2      100.000           3556.657
                                   ↑ each BED line measured independently (positions 1–200 for katG1, 150–400 for katG2)
```

The locus coverage report combines all reads from katG1 and katG2 into a single measurement across positions 1–400. The target report calculates coverage for each BED line separately.

### Resolving overlapping regions

In the example above, katG1 and katG2 share an overlapping region (positions 150–200). Without overlap resolution, a read spanning that region gets counted in **both** targets — inflating the locus coverage when they're aggregated together.

With `--resolve_overlapping_regions` enabled, tbp-parser identifies reads that **only** appear in non-overlapping portions of each target and uses those as a whitelist:

```
katG1: ├─────────────────┤
        1               200
              |~overlap~|
katG2:        ├───────────────────┤
              150               400

Non-overlapping katG1:   1 ─── 149  (unique to katG1)
Overlap zone:          150 ─── 200  (shared)
Non-overlapping katG2: 201 ─── 400  (unique to katG2)
```

For each target, tbp-parser will:

1. Finds reads in the **non-overlapping** portion of that target
2. Uses those read names as a whitelist
3. Reinterprets all positions (including the overlap zone), keeping **only** the whitelisted reads

This prevents double-counting when the targets are aggregated into locus coverage. The impact is visible in both reports:

```
WITHOUT --resolve_overlapping_regions:

locus_coverage_report.csv:
sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv1908c    katG       100.000           3471.463  ← inflated by double-counted reads

target_coverage_report.csv:
sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv1908c    katG1      100.000           2202.284  ← inflated by double-counted reads
sample01     Rv1908c    katG2      100.000           3556.657  ← inflated by double-counted reads
```

```
WITH --resolve_overlapping_regions:

locus_coverage_report.csv:
sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv1908c    katG       100.000            2683.554  ← only whitelisted reads

target_coverage_report.csv:
sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv1908c    katG1      100.000            2068.351  ← only whitelisted reads
sample01     Rv1908c    katG2      100.000            2766.803  ← only whitelisted reads
```

Note that with overlap resolution, coverage numbers may decrease in both reports. Reads that appear entirely within the overlapping region (no part spills into the non-overlapping regions) are excluded entirely since the target region origin cannot be conclusively determined. This is a conservative and more accurate approach for handling tNGS data with overlapping regions.

!!! warning "Trade-off: short reads and large overlaps"
    Overlap resolution works by whitelisting reads that appear in non-overlapping regions. If a read is short enough to exist **solely** within the overlap zone (i.e., it never extends into either target's unique region), it will be dropped from both targets — there is no way to determine which target it came from. This means that the larger the overlap relative to your read length, the more reads may be excluded. For BED files with large overlapping regions and short reads, this can noticeably reduce coverage numbers. Consider this trade-off when deciding whether to enable `--resolve_overlapping_regions`.

## Quick reference

| Scenario | Reports generated | Locus report | Target report |
|----------|-------------------|-------------|---------------|
| 1 BED line per gene | Locus only | 1 row per gene | Not generated |
| Multiple BED lines, same locus | Both | Aggregated into 1 row | 1 row per BED line |
| Overlapping BED lines + flag | Both | Aggregated (deduplicated) | 1 row per BED line (deduplicated) |
| Overlapping BED lines, no flag | Both | Aggregated (may double-count) | 1 row per BED line (may double-count) |

### Deletion warnings in coverage reports

When a variant contains a deletion in the ORF that passes positional QC and overlaps with a coverage region, it is flagged as a "valid deletion" and reported in the `qc_warning` column. These warnings explain why a gene's coverage may be lower than expected — the deleted region genuinely has no reads, which is biologically correct rather than a sequencing failure.

A deletion is considered valid if:

1. It is a deletion within the ORF (open reading frame)
2. It passes positional QC (sufficient depth, frequency, and read support)
3. Its genomic coordinates overlap with the coverage region

#### How `--err_coverage_bed` affects deletion warnings

The `qc_warning` column reports deletions differently depending on whether ERR coverage is active:

- **Without `--err_coverage_bed`**: All valid deletions in the full target/locus region are reported
- **With `--err_coverage_bed`**: Only valid deletions that fall **within the ERR region** are reported

A deletion that exists in the full target region but falls outside the ERR region will **not** appear in `qc_warning` when ERR coverage is active. The warning text also indicates which region was used:

```
Without ERR:  "Contains valid deletion(s) in full locus region: ..."
With ERR:     "Contains valid deletion(s) in ERR locus region: ..."
```

#### Example

Consider rpoB split across two target regions, with three deletions:

- `c.727_728delGTinsAC` — falls within the full target but **outside** the ERR region
- `c.1291_1292delAGinsCC` — falls within both the full target and the ERR region
- `c.1327_1329delTTGinsCCC` — falls within both the full target and the ERR region

```
Without --err_coverage_bed:

locus_coverage_report.csv:
sample_name  locus_tag  gene_name  ...  qc_warning
sample01     Rv0667     rpoB       ...  Contains valid deletion(s) in full locus region: c.727_728delGTinsAC; c.1291_1292delAGinsCC; c.1327_1329delTTGinsCCC

target_coverage_report.csv:
sample_name  locus_tag  gene_name  ...  qc_warning
sample01     Rv0667     rpoB_1     ...  Contains valid deletion(s) in full target region: c.727_728delGTinsAC
sample01     Rv0667     rpoB_2     ...  Contains valid deletion(s) in full target region: c.1291_1292delAGinsCC; c.1327_1329delTTGinsCCC
```

```
With --err_coverage_bed:

locus_coverage_report.csv:
sample_name  locus_tag  gene_name  ...  qc_warning
sample01     Rv0667     rpoB       ...  Contains valid deletion(s) in ERR locus region: c.1291_1292delAGinsCC; c.1327_1329delTTGinsCCC

target_coverage_report.csv:
sample_name  locus_tag  gene_name  ...  qc_warning
sample01     Rv0667     rpoB_1     ...                              ← no warning (deletion outside ERR)
sample01     Rv0667     rpoB_2     ...  Contains valid deletion(s) in ERR target region: c.1291_1292delAGinsCC; c.1327_1329delTTGinsCCC
```

Note that `c.727_728delGTinsAC` disappears from the warnings when ERR coverage is active because it falls outside the ERR region.

#### Why deletions matter for QC

Valid deletions also affect locus QC decisions. When a locus has low breadth of coverage (below `--min_percent_coverage`), a valid deletion in that locus can prevent the variant from being marked as "Insufficient Coverage" — the low coverage is explained by a real biological event rather than a sequencing failure. See the [interpretation documentation](../algorithm/interpretation.md) for more details.

**Note:** The `--use_err_as_brr` flag does **not** affect the coverage report. That flag controls whether ERR regions are used in place of full regions for locus QC decisions (e.g., determining "Insufficient Coverage"). The coverage report's `qc_warning` column is controlled solely by whether `--err_coverage_bed` is provided — if it is, only deletions within the ERR region are reported regardless of the `--use_err_as_brr` setting.

## Customizing column names

To overwrite any of the output column names or text in the coverage reports, please use the following format in a [configuration file](../inputs.md#configuration-file) or use the command-line parameter `--find_and_replace`:

```yaml
FIND_AND_REPLACE:
  "gene_name": "My_Gene_Column"
  "percent_coverage": "My_Percent_Coverage_Column"
  ...
```

Please note that this will rename every instance of that text in **all** output reports (every instance of "gene_name" will be renamed to "My_Gene_Column" in all output files, etc.).
