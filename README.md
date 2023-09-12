# tbp-parser

This repository contains the tbp-parser tool which parses the JSON output of [Jody Phelan's TBProfiler tool](https://github.com/jodyphelan/TBProfiler). Available as a download-able Python package and as a Docker image, tbp-parser converts the output of TBProfiler into four files:

1. A _Laboratorian_ report, which contains information regarding each mutation and its associated drug resistance profile in a CSV file. This file also contains two interpretation fields -- "Looker" and "MDL" which are generated using the CDC's expert rules for interpreting the severity of potential drug resistance mutations.
2. A _LIMS_ report, formatted specifically for STAR LIMS. This CSV report summarizes the highest severity mutations for each antimicrobial and lists the relevant mutations for each gene.
3. A _Looker_ report that contains condensed information suitable for generating a dashboard in Google's Looker studio.
4. A _coverage_ report that contains the percent coverage of each gene relative to the H37Rv reference genome in addition to any warnings, such as any deletions identified in the gene.

Please reach out to us at [theiagen@support.com](mailto:theiagen@support.com) if you would like any custom file formats and/or changes to these output files that suit your individual needs.

## Installation

### Docker

We recommend using the following Docker image to run tbp-parser:

```bash
docker pull us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:0.0.4
```

## Usage

```bash
usage: python /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]

Parses Jody Phelon's TBProfiler JSON output into three files:
- a Laboratorian report,
- a LIMS report
- a Looker report, and
- a coverage report

positional arguments:
  input_json
          the JSON file produced by TBProfiler
  input_bam
          the BAM file produced by TBProfiler

optional arguments:
  -h, --help
          show this help message and exit
  -v, --version
          show program's version number and exit
  -o, --output_prefix 
          the output file name prefix
          Do not include a space
          default="tb_parser"
  -d, --min_depth 
          the minimum depth of coverage to pass QC
          default=10
  -c, --coverage_threshold 
          the minimum percent coverage for a gene to pass QC
          default=100
  -s, --sequencing_method, -s 
          the sequencing method used to generate the data
          Enclose in quotes if includes a space
          default="Sequencing method not provided"
  -p, --operator 
          the operator who ran the sequencing
          Enclose in quotes if includes a space
          default="Operator not provided"
  --verbose
          increase output verbosity
  --debug
          increase output verbosity to debug;
          overwrites --verbose
```

Please note that the BAM file must have the accompanying BAI file in the same directory. It must also be named exactly the same as the BAM file but ending with a `.bai` suffix.

## Outputs

### Coverage report

Takes the naming convention of `<output_prefix>.percent_gene_coverage.csv`. This report contains the following fields:

- Gene
- Percent_Coverage
- Warning

### Laboratorian report

Takes the naming convention of `<output_prefix>.laboratorian_report.csv`. This report contains the following fields:

- sample_id: the sample name
- tbprofiler_gene_name: the gene name
- tbprofiler_locus_tag: the locus tag
- tbprofiler_variant_substitution_type: the variant substitution type (missense_variant, upstream_gene_variant...)
- tbprofiler_variant_substitution_nt: the nucleotide substitution (c.1349C>G)
- tbprofiler_variant_substitution_aa: the amino acid substitution (p.Ser450Trp)
- confidence: the tbprofiler annotation regarding resistance (Not assoc w R, Uncertain significance...)
- antimicrobial: the antimicrobial drug the mutation confers resistance to (streptomycin, rifampin...)
- looker_interpretation: the interpretation of resistance for the CDPH Looker report (R, R-interim, U, S, S-interim)
- mdl_interpretation: the MDL interpretation of resistance (R, S, U)
- depth: the depth of coverage at the mutation site (100)
- frequency: the frequency of mutation at the site (1)
- read_support: the number of reads supporting the mutation (100, depth*frequency)
- rationale: the rationale for resistance calling (WHO classification, Expert rule)
- warning: a column reserved for warnings such as low depth of coverage

### LIMS report

Takes the naming convention of `<output_prefix>.lims_report.csv`. This report contains the following fields:

- MDL sample accession numbers:  sample name
- M_DST_A01_ID - lineage
- The set of information in ANTIMICROBIAL_CODE_TO_GENES dictionary with target drug resistance information in layman's terms, and the mutations responsible for the predicted phenotype
- Date of analysis in YYYY-MM-DD HH:SS format
- Operator information

### Looker report

Takes the naming convention of `<output_prefix>.looker_report.csv`. This report contains the following fields:

- sample_id: the sample name
- for each antimicrobial, indication if resistant (R) or susceptible (S)
