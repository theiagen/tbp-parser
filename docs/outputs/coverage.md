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

Each line/entry in your BED file defines a **`BedRecord`** — a stretch of the genome where tbp-parser counts reads and calculates coverage. How you define these regions directly controls what gets measured and what appears in your coverage reports.

- The **locus** coverage report aggregates all `BedRecord`s that share the same locus tag and considers them as a single unit.
- The **target** coverage report keeps each `BedRecord` as its own separate unit.

!!! info "Which report gets generated?"
    If any `BedRecord`s share the same locus tag, both the **locus** and **target** coverage reports will be generated. Otherwise, only the **locus** coverage report will be generated.

!!! tip "QC uses locus-level coverage"
    All coverage-based QC determinations use the **locus-level** breadth of coverage.

    For example, if primer 1 for katG has poor coverage but primer 2 for katG has sufficient coverage, QC checks use the **overall** ratio of the aggregated coverage, not the individual primer regions.

---

### Example 1: _one target region in BED file_

When a BED file contains a single `BedRecord` (one **gene name** for one **locus tag**),

```
BED file:
Chrom    1        500    Rv0005       gyrB
```

only the locus coverage report is produced:

```
locus_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv0005     gyrB       100.000           542.310

  → percent_coverage = (num_positions ≥ min_depth) / num_positions
  → average_depth = (sum of reads at every position) / num_positions
```

---

### Example 2: _split regions_

A gene can be covered by multiple non-overlapping target regions (common in tNGS), resulting in multiple `BedRecord`s with the **same locus tag** but **different gene names** and **a gap** between them:

```
BED file:
Chrom    100      200    Rv0667       rpoB_1
Chrom    300      400    Rv0667       rpoB_2

├── rpoB_1 ──┤              ├── rpoB_2 ──┤
100         200            300         400
                    ↑
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
                                   ↑ each `BedRecord` measured independently
```

Note that the gap (positions 201–299) is **not measured at all** and no reads are counted there. The locus report aggregates only the positions defined by the `BedRecord`s (100–200 and 300–400), meaning locus-level breadth of coverage reflects only the regions defined in the BED file, not the entire gene.

---

### Example 3: _overlapping regions_

A gene can be covered by multiple overlapping target regions (common in tNGS), resulting in multiple `BedRecord`s with the **same locus tag** but **different gene names** that **overlap**:

```
BED file:
Chrom    1        200    Rv1908c      katG1
Chrom    150      400    Rv1908c      katG2

├──── katG1 ────┤
1               200
            ├────── katG2 ──────┤
            150                 400
```

Because there are more targets (2) than loci (1), **both** reports are generated:

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
                                   ↑ each `BedRecord` measured independently (positions 1–200 for katG1, 150–400 for katG2)
```

The locus coverage report combines all reads from katG1 and katG2 into a single measurement across positions 1–400. The target report calculates coverage for each `BedRecord` separately.

---

### Example 4: _resolving overlapping regions_

In [Example 3 above](#example-3-overlapping-regions), katG1 and katG2 share an overlapping region (positions 150–200). Without overlap resolution, a read spanning that region gets counted in **both** `BedRecord`s — inflating the locus coverage when they're aggregated together.

With `--resolve_overlapping_regions` enabled, tbp-parser identifies reads that **only** appear in _non_-overlapping portions of each `BedRecord` and uses those as a whitelist:

```
katG1: ├──────────────────┤
        1               200
              ├ overlap ┤
katG2:        ├───────────────────┤
              150               400

Non-overlapping katG1:   1 ─── 149  (unique to katG1)
Overlap zone:          150 ─── 200  (shared)
Non-overlapping katG2: 201 ─── 400  (unique to katG2)
```

For each `BedRecord`, tbp-parser will:

1. Find all reads that appear in the **non-overlapping** portion of that `BedRecord`
2. Use those read names as a whitelist
3. Reanalyze all positions (including the overlap zone), keeping **only** the whitelisted reads

This prevents double-counting when the `BedRecord`s are aggregated together in the locus coverage report. The impact is visible in both reports:

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

!!! warning "Trade-off: short reads and large overlaps"
    Note that overlap resolution may reduce coverage in both reports. Reads that fall entirely within an overlapping region, without extending into any unique region, are excluded because they cannot be confidently assigned to a specific `BedRecord`. This means that larger overlaps relative to the read length result in more reads being excluded. This is a conservative and more accurate approach for handling tNGS data with overlapping regions. Consider this trade-off when deciding whether to enable `--resolve_overlapping_regions`.

---

## QC warnings in coverage reports

A deletion is considered valid if:

1. It is a deletion within the ORF (open reading frame)
2. It passes positional QC (sufficient depth, frequency, and read support)
3. Its genomic coordinates overlap with the coverage region (either locus or target) being measured

Valid deletions are flagged and reported in the `qc_warning` column and may explain why a gene's coverage may be lower than expected.

!!! info "How `--err_coverage_bed` affects deletion warnings"
    The `qc_warning` column reports deletions differently depending on whether ERR coverage is active:

    - **Without `--err_coverage_bed`**: All valid deletions in the full target/locus region are reported
    - **With `--err_coverage_bed`**: Only valid deletions that fall **within the ERR region** are reported

    A deletion that exists in the full locus/target region but falls outside the ERR region will **not** appear in `qc_warning` when an `--err_coverage_bed` file is provided.

!!! info "The `--use_err_as_brr` flag does not affect the coverage report"
    The `--use_err_as_brr` flag controls whether ERR regions are used in place of full regions for locus QC decisions (e.g., determining "Insufficient Coverage"). The coverage report's `qc_warning` column is controlled solely by whether `--err_coverage_bed` is provided — if it is, only deletions within the ERR region are reported regardless of the `--use_err_as_brr` setting.

Deletions are shown only for the region they overlap. Consider rpoB split across two `BedRecord`s, with two deletions:

- `c.727_728delGTinsAC` — falls within rpoB_1
- `c.1291_1292delAGinsCC` — falls within rpoB_2

In the **locus** coverage report, all deletions across both `BedRecord`s are aggregated into a single warning for the locus:

```
locus_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv0667     rpoB       85.000            280.500        Contains valid deletion(s) in full locus region: c.727_728delGTinsAC; c.1291_1292delAGinsCC
```

In the **target** coverage report, each deletion only appears on the specific `BedRecord` it overlaps with:

```
target_coverage_report.csv:

sample_name  locus_tag  gene_name  percent_coverage  average_depth  qc_warning
sample01     Rv0667     rpoB_1     90.000            350.100        Contains valid deletion(s) in full target region: c.727_728delGTinsAC
sample01     Rv0667     rpoB_2     80.000            210.900        Contains valid deletion(s) in full target region: c.1291_1292delAGinsCC
```

The target coverage report lets you see exactly which region contains which deletion, while the locus coverage report gives the combined picture.

---

## Customizing column names

To overwrite any of the output column names or text in the coverage reports, please use the following format in a [configuration file](../inputs.md#configuration-file) or use the command-line parameter `--find_and_replace`:

```yaml
FIND_AND_REPLACE:
  "gene_name": "My_Gene_Column"
  "percent_coverage": "My_Percent_Coverage_Column"
  ...
```

Please note that this will rename every instance of that text in **all** output reports (every instance of "gene_name" will be renamed to "My_Gene_Column" in all output files, etc.).
