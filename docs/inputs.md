---
title: Command-line Arguments
---

This page documents the arguments of the main `parse` subcommand, which produces the four output reports.

## Required Arguments

`tbp-parser` is designed to run immediately after [Jody Phelan’s TBProfiler tool](https://github.com/jodyphelan/TBProfiler). Five inputs are required:

| Parameter  | Description | Purpose |
| :--------- | :---------- | :------ |
| `--input_json` | The path to the results JSON file that was produced by `TBProfiler` v6+ | Contains information about the mutations detected in the sample: quality, type, and any antimicrobial resistance information. |
| `--input_bam`  | The path to the BAM file that was produced by `TBProfiler` v6+ | Contains the alignment information for the sample; needed for determining sequencing quality for quality control. Please note that the BAM file should have the accompanying BAI file in the same directory; if it is missing, `tbp-parser` will generate one, which can take a while. |
| `--coverage_bed` | The path to a BED file containing the genes of interest, their locus tags, and their regions ([see below](#coverage-bed-file)) | Defines the regions that breadth of coverage and average depth are calculated over. |
| `--gene_database_yml` | The path to a gene database YAML file ([see below](#gene-database-file)) | Represents the TBProfiler database the variants were called against: which genes exist, which drugs they are associated with, and each gene's tier and promoter region. |
| `--lims_report_format_yml` | The path to a LIMS report format YAML file ([see below](#lims-report-format-yaml-file)) | Defines the columns of the LIMS report. |

!!! dna "No default input files are provided"
    As of v4.0.0, `tbp-parser` no longer ships with default versions of these files, because the correct contents depend on which TBProfiler database your sample was called against. Generate the gene database and LIMS report format from that database's `genes.bed` file with [`tbp-parser build_gene_db` and `tbp-parser build_lims_fmt`](subcommands.md), then reuse them for every sample called against that database.

### Coverage BED File

The Coverage BED file is the **tab-delimited** [BED](https://grch37.ensembl.org/info/website/upload/bed.html) file that contains gene regions of interest and their associated antimicrobials. This file is used for quality control calculations. The file should be formatted like the [genes.bed](subcommands.md#where-to-get-a-genesbed) file in TBProfiler, with the following columns in this order:

1. `chrom`: the chromosome or contig name which must match the chromosome name in the BAM file (e.g. "Chromosome")
2. `start`: the start position of the gene (e.g. 1)
3. `end`: the end position of the gene (e.g. 1524)
4. `locus_name`: the locus name of the gene (e.g. "Rv0001")
5. `gene_name`: the gene name (e.g. "dnaA")
6. `drugs`: the drugs associated with that gene, separated by commas (e.g. "isoniazid,rifampicin")

??? warning "The sixth `drugs` column"
    - The `parse` (`--coverage_bed`) subcommand reads only the first five columns; `drugs` and anything after it are ignored, because the gene database is the truth set for gene-drug associations, not the coverage BED file.
    - The `build_gene_db` (`--db_bed`) subcommand is the opposite: it builds the database *from* this column, and **a gene with no drugs listed is left out of the database entirely**.
    - Populating this column in every BED file means the same file can serve both commands.

For example, the following is a valid BED file:

```text
Chromosome	1	    1524	Rv0001	dnaA	isoniazid
Chromosome	4933	7267	Rv0005	gyrB	levofloxacin,moxifloxacin
```

Please note that this file *does not* have a header line. [See the Subcommands page for where to obtain a `genes.bed` file](subcommands.md#where-to-get-a-genesbed) for the TBProfiler database you ran your sample against.

This is also the same format used for the optional `--err_coverage_bed` file, which is an optional input parameter primarily for tNGS analysis ([see below](#tngs-specific-arguments)).

---

### Gene Database File

The *gene database* file contains a dictionary of the following information for each gene:

1. `locus_tag`: the locus tag of the gene (e.g. `Rv0005`)
2. `gene_name`: the gene name (e.g. `gyrB`)
3. `tier`: the tier of the gene (e.g. `Tier 1`)
4. `promoter_region`: the WHO-specified proximal promoter region (e.g. `[-108, -1]`)
5. `drugs`: the antimicrobials associated with this gene (e.g. `[levofloxacin, moxifloxacin]`)
6. `aliases` *(optional)*: alternate locus tags the gene may be referred to by (e.g. `[Rvnr01, MTB000019]` for `rrs`)

This file represents the TBProfiler database your variants were called against, and is [generated with `tbp-parser build_gene_db`](subcommands.md#build_gene_db) from that database's `genes.bed` file. If you would like to include a different gene, or modify the content of existing entries, you can do so by using the following format:

```yaml
# do not modify unbracketed text
# text within angle brackets should be replaced with the appropriate information for the gene of interest
<locus_tag_of_gene>:
  locus_tag: <locus_tag_of_gene>
  gene_name: <gene_name>
  tier: <tier_of_gene>
  promoter_region: [<WHO-specified_proximal_promoter_regions_start>, <WHO-specified_proximal_promoter_regions_end>]
  drugs: [<drug_1>, <drug_2>, ...]
<locus_tag_of_gene2>:
  ...
...
```

If information for your gene of interest is not available, please use the following values as placeholders:

- for `tier`, use `NA`
- for `promoter_region`, use `[]`

A gene may also have two promoter windows, in which case `promoter_region` is a nested list (e.g. `[[-51, -1], [-503, -323]]`).

For example, the following are valid entries in the gene database file:

```yaml
Rv0001:
    locus_tag: Rv0001
    gene_name: dnaA
    tier: NA
    promoter_region: [-314, -1]
    drugs: [isoniazid]
Rv0676c:
    locus_tag: Rv0676c
    gene_name: mmpL5
    tier: Tier 1
    promoter_region: []
    drugs: [bedaquiline, clofazimine]
EBG00000313325:
    locus_tag: EBG00000313325
    gene_name: rrs
    tier: Tier 1
    promoter_region: [-151, -1]
    drugs: [amikacin, capreomycin, kanamycin, streptomycin]
    aliases: [Rvnr01, MTB000019]
```

See [`build_gene_db`](subcommands.md#build_gene_db) for more information on how to determine the correct values for those fields.

---

### LIMS Report Format YAML File

Different LIMS systems may require different column formatting for easy import. The LIMS report format YAML file allows users to specify the output column names for the LIMS report output. This file is required; the [`build_lims_fmt`](subcommands.md#build_lims_fmt) subcommand will generate one containing every gene-drug combination in your gene database, which you can then edit.

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
[Please see the LIMS report section for more information on this input file, the report, its purpose, and additional customization options.](outputs/lims.md)

---

## Validation Arguments

This option controls whether `tbp-parser` checks that your input files agree with the gene database before doing any work.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| `--skip_input_validation` | Skip validation that all genes and gene/drug pairs referenced in the input files are present in the gene database | false |

Input files that cover fewer genes or drugs than the gene database are fine; anything the gene database does not know about is not, because it cannot be interpreted or reported on. Four checks are performed:

1. Every gene/drug pair in the --input_json exists in the gene database
2. Every gene/drug pair in the --lims_report_format_yml exists in the gene database
3. Every gene in the --lims_report_format_yml has a region defined in the --coverage_bed file
4. Every locus tag in the --coverage_bed exists in the gene database

Passing `--skip_input_validation` disables **only** these four checks. File accessibility, BED column count, BAM index, threshold range, boundary format, duplicate BED region, and [ERR region containment](#tngs-specific-arguments) checks all still run.

---

## Quality Control Arguments

These options determine the thresholds for quality control.

| Long Version           | Description | Default Value |
| :--------------------- | :---------- | :------------ |
| `--min_depth` | The minimum depth of coverage required for a site to pass QC | 10 |
| `--min_percent_coverage` | The minimum fraction of a region that has depth meeting or exceeding the threshold set by `min_depth` (used for a gene/locus to pass QC; 1.0 -> 100%) | 1.0 |
| `--min_read_support` | The minimum read support for a mutation to pass QC | 10 |
| `--min_frequency` | The minimum frequency for a mutation to pass QC (0.1 -> 10%) | 0.1 |
| `--min_percent_loci_covered` | The minimum fraction of loci/genes in the LIMS report that must pass coverage QC for the sample to be identified as MTBC (0.7 -> 70%) | 0.7 |

---

## tNGS-specific Arguments

These options are primarily used for tNGS data.

| <div style="width:250px">Name</div> | Description | Default Value |
| :--- | :---------- | :------------ |
| `--tngs` | Indicates that the input data was generated using a tNGS protocol. Turns on tNGS-specific features | false |
| `--err_coverage_bed` | the BED file containing the "essential for resistance regions." This file indicates to tbp-parser that these regions should also have breadth of coverage and average depth calculations performed; this file should be formatted like the genes.bed file in TBProfiler and the [coverage BED described above](#coverage-bed-file) | |
| `--use_err_for_qc` | if an ERR BED file is provided, use the ERR coverage regions in place of the typical coverage regions for all QC determinations. This option is **experimental**.<br>Note: This will influence how variants are interpretated and how deletions are reported because the QC thresholds for breadth of coverage and average depth will be based on the coverage found within the ERR regions. | false |
| `--resolve_overlapping_regions` | Resolve overlapping BED regions to avoid double-counting reads across overlapping targets. Recommended for tNGS data with overlapping amplicon regions. See [Handling overlapping primer regions](./algorithm/technical.md#2-coverage-calculations) and the section on [the coverage report](./outputs/coverage.md) for more details | false |
| `--tngs_frequency_boundaries` | the frequency boundaries (comma-delimited; `lower_f,upper_f`) for tNGS QC reporting, used in conjunction with `--tngs_read_support_boundaries` | 0.1,0.1 |
| `--tngs_read_support_boundaries` | the read support boundaries (comma-delimited; `lower_rs,upper_rs`) for tNGS QC reporting, used in conjunction with `--tngs_frequency_boundaries` | 10,10 |

??? warning "ERR regions must fall inside their coverage regions"
    Each region in `--err_coverage_bed` must be entirely contained within the corresponding `--coverage_bed` region for the same gene/locus. If it is not, `tbp-parser` stops with an error naming the offending locus and gene, for example:

    ```text
    ERR coords for target `Rv0006` (gyrA) [(7000, 7500)] fall outside target coords [(7068, 9818)]
    ```

    A gene present in `--coverage_bed` but **absent** from `--err_coverage_bed` is not an error; it reports `N/A` in the [coverage report](outputs/coverage.md).

---

## Text Arguments

These options are used verbatim in the reports, or are used to name the output files.

| <div style="width:200px">Long Version</div> | Description | Default Value |
| :--- | :---------- | :------------ |
| `--sequencing_method` | The sequencing method used to generate the data; used in the LIMS & Looker reports. Enclose in quotes if including a space | "Sequencing method not provided" |
| `--operator` | The operator who ran the analysis; used in the LIMS & Looker reports. Enclose in quotes if including a space | "Operator not provided" |
| `--output_prefix` | The prefix to use for the output files; a trailing slash is treated as a directory, which is created if it does not exist. Do not include any spaces | "tbp_parser" |
| `--find_and_replace` | A JSON string that can be used to specify any text in the output files that should be find-and-replaced with other text. The keys will be the text to find, and the values will be the text to replace it with. This is useful for labs that want to customize the text in their reports (e.g. renaming drugs or genes or output columns).<br>For example, `'{"rifampicin": "rifampin", "fbiD": "Rv2983", "mmpR5": "Rv0678", "p.0?": ""}'` | '{}' |

---

## Logging Arguments

These options change the verbosity of the `stderr` log. The `parse` subcommand additionally writes the same log to `{output_prefix}.log`.

| Name | Description | Default Value |
| :--- | :---------- | :------------ |
| `--debug` | The highest level of output verbosity detailing every step of the analysis and logic implemented | false |

---

## Configuration File

Instead of providing the input parameters on the command line, the ability to provide a configuration file in YAML format is available. This file (and any included fields) are **case-sensitive** and should be written in all caps.

The configuration file will accept input parameters from the [Validation Arguments](#validation-arguments), [Quality Control Arguments](#quality-control-arguments), [Text Arguments](#text-arguments), and the [tNGS-specific Arguments](#tngs-specific-arguments). Any flag that requires a file path as an input should be provided **separately** on the command line. This includes: [Required Arguments](#required-arguments), `--err_coverage_bed`, and [Logging Arguments](#logging-arguments). The configuration file can be provided using the `--config` argument. Input parameters should be indicated in all caps and should match the long version of the command-line arguments (e.g. `MIN_FREQUENCY` instead of `--min_frequency`; see example below).

```yaml
# I can overwrite any input parameters, like so.
# This makes it easy to rerun the same analysis on different
# samples without rewriting all of the parameters each time.
MIN_FREQUENCY: 0.1
MIN_PERCENT_LOCI_COVERED: 0.7
SKIP_INPUT_VALIDATION: false
TNGS: true
RESOLVE_OVERLAPPING_REGIONS: true
TNGS_FREQUENCY_BOUNDARIES:
- 0.1
- 0.95
TNGS_READ_SUPPORT_BOUNDARIES:
- 100
- 500

# I can also use the configuration file to customize output files.
# My laboratory reports "rifampicin" as "rifampin", so I want to
# rename that text in all of the output files. I also use Rv0678
# instead of mmpR5 and Rv2983 instead of fbiD; and I need to rename
# an output column in the LIMS report from "Sample Name" to "sample"
FIND_AND_REPLACE:
  rifampicin: "rifampin"
  fbiD: "Rv2983"
  mmpR5: "Rv0678"
  "Sample Name": "sample"
```

---