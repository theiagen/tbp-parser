# tbp-parser

**CAUTION:** The information produced by this program should **not** be used for clinical reporting unless and until extensive validation has occured in your laboratory on a stable version. Otherwise, the outputs of tbp-parser are for research use only.

## Overview

This repository contains the tbp-parser tool which parses the JSON output of [Jody Phelan's TB-Profiler tool](https://github.com/jodyphelan/TBProfiler). Available as a downloadable Python package and as a Docker image, tbp-parser converts the output of TB-Profiler into four files:

1. A _Laboratorian_ report, which contains information regarding each mutation and its associated drug resistance profile in a CSV file. This file also contains two interpretation fields -- "Looker" and "MDL" which are generated using the CDC's expert rules for interpreting the severity of potential drug resistance mutations.
2. A _LIMS_ report, formatted specifically for STAR LIMS. This CSV report summarizes the highest severity mutations for each antimicrobial and lists the relevant mutations for each gene.
3. A _Looker_ report that contains condensed information suitable for generating a dashboard in Google's Looker studio.
4. A _coverage_ report that contains the percent coverage of each gene relative to the H37Rv reference genome in addition to any warnings, such as any deletions identified in the gene.

Please reach out to us at [support@theiagen.com](mailto:support@theiagen.com) if you would like any custom file formats and/or changes to these output files that suit your individual needs.

See also [this page](https://theiagen.notion.site/tbp-parser-b02bef0cbc814b129875d861698c80a2?pvs=4) for additional documentation on how tbp-parser works and where the expert rules were derived.

## Brief Description of Versions

You may notice there are many releases; tbp-parser is in active development and each release is "use at your own risk." We highly recommend upgrading to the latest release as they include important bug fixes. In order to help track the different changes, we have included a brief description of each release:

- **v1.2.x _& below_** - the initial developmental stages of tbp-parser for WGS data
- **v1.3.x** - the addition of tNGS data parsing and includes some updates applicable to WGS parsing
- **v1.4.x** - reworks how QC is performed (changes in order of operations)
  - **v1.4.3+** - changes how tNGS lineage determination is performed
  - **v1.4.4+** - changes how nonsynonymous mutations are interpretted; major interpretation differences between earlier versions
- **v1.5.x** - major changes to code in that it expects results from TB-Profiler v6.2.0+; no longer backwards compatible and really should've been a v2 release but it's too late now
  - code changes for v1.5x are available on the `who-v2` branch

Again, please use tbp-parser at your own risk and be sure to perform extensive validation before using this tool in a clinical setting.

## Installation

### Docker

We highly recommend using the following Docker image to run tbp-parser:

```markdown
docker pull us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:1.4.4.2
```

The entrypoint for this Docker image is the tbp-parser help message. To run this container interactively, use the following command:

```markdown
docker run -it --entrypoint=/bin/bash us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:1.4.4.2
# Once inside the container interactively, you can run the tbp-parser tool
python3 /tbp-parser/tbp_parser/tbp_parser.py -v
# v1.4.4.2
```

### Locally with Python

tbp-parser is not yet available with `pip` or `conda`. To run tbp-parser in your local command-line environment, install the following dependencies:

- python3
- pandas >= 1.4.2
- importlib_resources
- samtools

After installation, download and extract the latest release of tbp-parser and run the script with Python.

## Usage

The following is the help message printed by tbp-parser.

```markdown
usage: python3 /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]

Parses Jody Phelon's TB-Profiler JSON output into three files:
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

quality control arguments:
  options that determine what passes QC

  -d, --min_depth
          the minimum depth of coverage for a site to pass QC
          default=10
  -c, --min_percent_coverage
          the minimum percentage of a region that has depth above the threshold set by min_depth
            (used for a gene/locus to pass QC)
          default=100
  -s, --min_read_support
          the minimum read support for a mutation to pass QC
          default=10
  -f, --min_frequency
          the minimum frequency for a mutation to pass QC (0.1 -> 10%)
          default=0.1
  -r, --coverage_regions
          the BED file containing the regions to calculate percent coverage for
          default=data/tbdb-modified-regions.bed

text arguments:
  arguments that are used verbatim in the reports or to name the output files

  -m, --sequencing_method
          the sequencing method used to generate the data; used in the LIMS & Looker reports
          ** Enclose in quotes if includes a space
          default="Sequencing method not provided"
  -p, --operator
          the operator who ran the sequencing; used in the LIMS & Looker reports
          ** Enclose in quotes if includes a space
          default="Operator not provided"
  -o, --output_prefix
          the output file name prefix
          ** Do not include any spaces

tNGS-specific arguments:
  options that are primarily used for tNGS data
  (all frequency arguments are compatible with WGS data)

  --tngs
          indicates that the input data was generated using Deeplex + CDPH modified protocol
          Turns on tNGS-specific global parameters
  --tngs_expert_regions
          the BED file containing the regions to calculate coverage for expert rule regions
            (used to determine coverage quality in the regions where resistance-conferring
            mutations are found, or where a CDC expert rule is applied; not for QC)
          default=data/tngs-expert-rule-regions.bed
  --rrs_frequency
          the minimum frequency for an rrs mutation to pass QC
            (rrs has several problematic sites in the Deeplex tNGS assay)
          default=0.1
  --rrl_frequency
          the minimum frequency for an rrl mutation to pass QC
            (rrl has several problematic sites in the Deeplex tNGS assay)
          default=0.1
  --rpob449_frequency
          the minimum frequency for an rpoB mutation at protein position 449 to pass QC
            (this is a problematic site in the Deeplex tNGS assay)
          default=0.1
  --etha237_frequency
          the minimum frequency for an ethA mutation at protein position 237 to pass QC
            (this is a problematic site in the Deeplex tNGS assay)
          default=0.1

logging arguments:
  options that change the verbosity of the stdout log

  --verbose
          increase output verbosity
  --debug
          increase output verbosity to debug; overwrites --verbose

Please contact support@theiagen.com with any questions
```

Please note that the BAM file must have the accompanying BAI file in the same directory. It must also be named exactly the same as the BAM file but ending with a `.bai` suffix.

### Example

This shows how the script can be run if used inside the Docker container provided above.

```markdown
python3 /tbp-parser/tbp_parser/tbp_parser.py \
    /path/to/data/tbprofiler_output.json \
    /path/to/data/tbprofiler_output.bam \
    -o "example" \
    --min_depth 12 \
    --min_frequency 0.9 \
    --sequencing_method "Illumina NextSeq" \
    --operator "John Doe" 
```

## Outputs

### Coverage report

Takes the naming convention of `<output_prefix>.percent_gene_coverage.csv`. This report contains the following fields:

- **Gene** - the name of the gene or locus
- **Percent_Coverage** - the percent of the gene (positions determined by the H37Rv genome) that is covered at a depth greater than the `--min_depth` value
- **Warning** - indicates if any deletions were identified in the gene which may contribute to lower than expected coverage

If the `--tngs` flag is used, the report contains the following fields:

- **Gene** - the name of the gene or locus
- **Coverage_Breadth_reportableQC_region** - the percent of the gene (positions determined by the regions covered by the tNGS Deeplex + CDPH assay primers that are considered reportable by CDPH) that is covered at a depth greater than the `--min_depth` value
- **QC_Warning** - indicates if any deletions were identified in the gene which may contribute to lower than expected coverage
- **Coverage_Breadth_R_expert-rule_region** - the percent of the regions (positions that could contain any resistance-conferring mutations or require expert-rule application) that is covered at a depth greater than the `--min_depth` value

Coverage regions are determined with either the default "../data/tbdb-modified-regions.bed" (collected on Sep 1, 2023 from the TBProfiler repository, or if `--tngs`, "../data/tngs-reportable-regions.bed". 

The R-expert rule region is determined only if `--tngs` is indicated and uses the ranges in "../data/tngs-expert-rule-regions.bed". 

User-provided coverage regions always take precedence over default values

### Laboratorian report

Takes the naming convention of `<output_prefix>.laboratorian_report.csv`. This report contains the following fields:

- **sample_id**: the sample name
- **tbprofiler_gene_name**: the gene name
- **tbprofiler_locus_tag**: the locus tag
- **tbprofiler_variant_substitution_type**: the variant substitution type (missense_variant, upstream_gene_variant, ...)
- **tbprofiler_variant_substitution_nt**: the nucleotide substitution (e.g., c.1349C>G)
- **tbprofiler_variant_substitution_aa**: the amino acid substitution (e.g., p.Ser450Trp)
- **confidence**: the tbprofiler annotation regarding resistance (Not assoc w R, Uncertain significance, ...)
- **antimicrobial**: the antimicrobial drug the mutation confers resistance to (streptomycin, rifampin, ...)
- **looker_interpretation**: the interpretation of resistance for the CDPH Looker report (R, R-interim, U, S, or S-interim)
- **mdl_interpretation**: the MDL interpretation of resistance (R, S, or U)
- **depth**: the depth of coverage at the mutation site (e.g., 100)
- **frequency**: the frequency of mutation at the site (e.g., 1)
- **read_support**: the number of reads supporting the mutation (e.g., 100, calculated by `depth*frequency`)
- **rationale**: the rationale for resistance calling (e.g., WHO classification or expert rule)
- **warning**: a column reserved for warnings such as low depth of coverage or poor read support

### LIMS report

Takes the naming convention of `<output_prefix>.lims_report.csv`. This report contains the following fields:

- **"MDL sample accession numbers"** - sample name
- **M_DST_A01_ID** - the lineage of the sample in plain English
- The set of information in ANTIMICROBIAL_CODE_TO_GENES_WGS dictionary (or, if `--tngs`, ANTIMICROBIAL_CODE_TO_GENES_tNGS) with target drug resistance information in layman's terms, and the mutations responsible for the predicted phenotype
- **"Analysis date"** - the date of analysis in YYYY-MM-DD HH:SS format
- **Operator** - operator information, set by the `--operator` option
- **M_DST_O01_LINEAGE** - the `main_lin` of the sample as reported by TBProfiler

### Looker report

Takes the naming convention of `<output_prefix>.looker_report.csv`. This report contains the following fields:

- **sample_id**: the sample name
- **output_seq_method_type** - the sequencing method used to generate the data, set by the `--sequencing_method` option
- for each antimicrobial, indication if resistant (R), uncertain (U), or susceptible (S)
- **lineage**: the `main_lin` of the sample as reported by TBProfiler
- **ID** - the lineage of the sample in plain English (the same as `M_DST_A01_ID` in the LIMS report)
- **analysis_date** - the date of analysis in YYYY-MM-DD HH:SS format
- **operator** - operator information, set by the `--operator` option

## The interpretation logic document and where you can find the corresponding code

The interpretation logic document is long and complicated, and the implementation of these rules is scattered throughout the code. In order to make it easier for a newcomer to find things, I've listed the rules and where they are implemented in the code:

### Rule 1

- 1.1 - (WHO annotation) - `complete_row()` in `Row` (lines 245 - 274)
- 1.2 - (No WHO annotation; CDC expert rules) - `apply_expert_rules()` in `Variant` (lines 127 - 167)
- 1.3 - (mmpR5/Rv0678, mmpL5, and mmpS5 reporting) - `iterate_section()` in `Laboratorian` (lines 59 - 74), `extract_alternate_consequences()` in `Variant` (lines 37 - 58)

### Rule 2

- 2.1 - (WHO annotation) - `complete_row()` in `Row` (lines 258 - 262)
- 2.2 - (No WHO annotation; CDC expert rules) - `apply_expert_rules()` in `Variant` (lines 173 - 194, 215 - 219, 229 - 240)
  - 2.2.1 - (loss-of-function expert rule) - `apply_expert_rules()` in `Variant` (lines 176 - 194)
  - 2.2.2.1 - (rpoB RRDR expert rule) - `apply_expert_rules()` in `Variant` (lines 209 - 219)
  - 2.2.2.2 - (rpoB non-RRDR expert rule) - `apply_expert_rules()` in `Variant` (lines 229 - 235, 238 - 240)

### Rule 3

- 3.1 - (WHO annotation) - `complete_row()` in `Row` (lines 264 - 269)
- 3.2 - (No WHO annotation; CDC expert rules) - `apply_expert_rules()` in `Variant` (lines 190 - 240)
  - 3.2.1 - (rrs gene expert rule) - `apply_expert_rules()` in `Variant` (lines 245 - 253)
  - 3.2.2 - (GyrA QRDR expert rule) - `apply_expert_rules()` in `Variant` (lines 221 - 223)
  - 3.2.3 - (GyrB QRDR expert rule) - `apply_expert_rules()` in `Variant` (lines 225 - 227)
  - 3.2.4 - (All else expert rule) - `apply_expert_rules()` in `Variant` (lines 229 - 233, 236 - 240, 256 - 264)

### Rule 4

After here, it gets a little iffy because I haven't updated it in a while and it's exhausting try to track them down; the line numbers are more like general regions than exact coordinates at this point

- 4.1 - (WT gene has coverage) - `__init__()` in `Row` (lines 161 - 167, 182 - 188)
- 4.2 - (QC filtering) - `__init__()` in `Row` (lines 77 - 124, 133 - 135, 169 - 175, 190 - 197) and in the `create_laboratorian_report()` in `Laboratorian` (lines 131 - 157)
  - 4.2.1.1 - (Gene coverage QC pass) - `__init__()` in `Row` (lines 127 - 129)
  - 4.2.1.2 - (Gene coverage QC fail, but deletion) - `__init__()` in `Row` (lines 117 - 121, 133 - 135)
  - 4.2.1.3 - (Gene coverage position QC fail, no deletion) - `__init__()` in `Row` (lines 72 - 73, 128 - 142) and the `create_laboratorian_report()` in `Laboratorian` (lines 118 - 153)
    - 4.2.1.3.1 - (WT) - `__init__()` in `Row` (lines 128 - 142)
    - 4.2.1.3.2 - (No R mutations detected) - `create_laboratorian_report()` in `Laboratorian` (lines 118 - 153)
    - 4.2.1.3.3 - (R mutation detected) - `create_laboratorian_report()` in `Laboratorian` (lines 118 - 132)
  - 4.2.2 - (Variant position QC fail) - `__init__()` in `Row` (lines 83 - 92)
    - 4.2.2.1 (Variant position QC fail treated as WT) - `apply_lims_rules()` in `LIMS` (lines 180-186)
- 4.3 - (Additional notes) - `apply_expert_rules()` in `Variant` (lines 131 - 208) and `extract_annotation()` in `Variant` (lines 100 - 112)
  - 4.3.1 - (Indels spanning into regions of interest) - `apply_expert_rules()` in `Variant` (lines 131 - 134, 153-155, 190-208)
  - 4.3.2 - (When an SME expert is needed) - Is NOT implemented in the code as requires human intervention
  - 4.3.3 - (Report all antimicrobial drugs associated with a gene) - `extract_annotation()` in `Variant` (lines 100 - 112)

### Rule 5

- 5.1 - (Laboratorian language to LIMS language) - `apply_lims_rules()` in `LIMS` (line 180 - 182, 203 - 253)
- 5.2 - (Drug interpretation and lineage in LIMS report) - `apply_lims_rules()` in `LIMS` (line 108, 312) and `get_id()` in `LIMS` (line 49-50)
  - 5.2.1 - (Only report LIMS genes in the LIMS report) - `apply_lims_rules()` in `LIMS` (line 108)
  - 5.2.2 - (Report the "main_lin" field from TBProfiler) - `get_id()` in `LIMS` (line 49-50) and the `create_lims_rules()` in `LIMS` (line 312)
- 5.3 - (Individual mutation output format in the LIMS report) - `apply_lims_rules()` in `LIMS` (lines 132 - 196)
  - 5.3.1 - (Only report "R" or "U" mutations except rpoB RRDR) - `apply_lims_rules()` in `LIMS` (lines 185 - 190)
  - 5.3.2 (Mutation output format) - `apply_lims_rules()` in `LIMS` (lines 132, 160 - 196)
  - 5.3.3 (Take highest read support mutation) - `apply_lims_rules()` in `LIMS` (lines 133 - 157)
- 5.4 - (ID formatting in the LIMS report) - `get_id()` in `LIMS` (lines 25 -74)
  - 5.4.1 - (QC check) - `get_id()` in `LIMS` (lines 33 - 38, 53)
    - 5.4.1.1 - 5.4.1.4 - (Lineage formatting) - `get_id()` in `LIMS` (lines 54 - 65)
  - 5.4.2 - (QC fail) - `get_id()` in `LIMS` (lines 67 - 68)
  - 5.4.3 - (Reporting requirements) - Is NOT implemented in the code as requires human intervention
  - 5.4.4 - (Resequencing requirements) - Is NOT implemented in the code as requires human intervention

### Rule 6

- 6.1.1 - (Report all genes in the Laboratorian report for the Looker report) - `create_looker_report()` in `Looker` (line 63)
- 6.1.2 - (ID field for the Looker report) - `create_looker_report()` in `Looker` (line 57)
- 6.1.3 - (Lineage field for the Looker report) - `create_looker_report()` in `Looker` (line 56)

## Description of Algorithm

For those morbidly interested in how tbp-parser works, you can read a step-by-step explanation of how the TBProfiler outputs are parsed below.

1. The `tbp_parser.py` script is run with the input JSON and BAM files. This script is the entry point for the tool, and begins by parsing the various command-line arguments. Upon completion of arguments, it creates an new `Parser` object with those arguments as input and runs the `.run()` method of that class.
2. The Parser class (`Parser.py`) is the controlling module for this package. It is initialized with the command line arguments and uses these inputs to set global variables, such as `MIN_DEPTH`, `COVERAGE_THRESHOLD`, etc. Its primary method is `run()`. This method also calls the `check_dependency_exists()` function, which determines if the prequesite `samtools` is installed. 
3. The first object created in the `run()` function is an instance of the `Coverage` class (`Coverage.py`) in order to make the initial coverage report. This class is initialized with the input BAM file. It then runs the `create_coverage_report()` method, which creates the coverage report and saves it to a CSV file.
    - The `calculate_coverage()` method iterates through the genes contained in the `../data/tbdb-modified-regions.bed` (or a different bed file provided by the user or as indicated by the `--tngs` flag) file (collected on Sep 1, 2023 from the TBProfiler repository). It runs `samtools depth` and then divides the sum of the positions above the `MIN_DEPTH` value by the length of the gene in question to calculate the percent coverage of each gene in the H37Rv reference genome. It then saves this information as `{gene: percent coverage}` into the global variable `COVERAGE_DICTIONARY`.
4. The section object created in the `run()` function is of the `Laboratorian` class (`Laboratorian.py`). This class is initialized with the input JSON file. The `create_laboratorian_report()` method is then called on the object, which begins by performing the `iterate_section()` method on the different variant sections belonging to its input JSON.
    - The `iterate_section()` takes two items: the variant section of the input JSON and the `row_list` which is a list of `Row` objects (explained later) to be included in the Laboratorian report. The method iterates through each variant in the section and creates a `Variant` object (`Variant.py`) for each one.
    - The `Variant` class represents a single variant reported by TBProfiler. This class is initialized by turning each item of the variant's dictionary into an attribute. The first method called on the object is `extract annotations()` which divides the annotation field of a `Variant` into individual annotations. Depending on if the object has an `annotation` attribute, each item in that field is turned into a `Row` object (`Row.py`). If not, every antimicrobial drug in the `gene_associated_drugs` attribute is transformed into a `Row`. This is because we preferentially examine the `annotation` field as it may contain a WHO resistance interpretation which takes priority over any expert rules applied in the `gene_associated_drugs` field.
    - The `Row` class represents a single row in the Laboratorian report; each attribute is a column in the report. Its initialization step is fairly complicated. As input, it can take the originating `Variant` object, a `who_confidence`, the associated antimicrobial drug, the originating gene (default _None_), the depth of the mutation (default 0), and the frequency of the mutation (default _None_).
        - If the originating `Variant` object is _not_ None, the `Row` object is initialized with the following attributes taken directly from the originating `Variant`:
            - `tbprofiler_gene_name`: the gene name (optionally provided during initialization)
            - `tbprofiler_locus_tag`: the locus tag
            - `tbprofiler_variant_substitution_type`: the type of mutation
            - `tbprofiler_variant_substitution_nt`: the mutation in nucleotide notation
            - `tbprofiler_variant_substitution_aa`: the mutation in amino acid notation
            - `confidence`: The WHO confidence ranking for the mutation (optionally provided during initialization)
            - `antimicrobial`: the antimicrobial drug the mutation/gene may be associated with (optionally provided during initialization)
            - `looker_interpretation`: the interpretation of resistance for the CDPH Looker report
            - `mdl_interpretation`: the MDL interpretation of resistance for the CDPH LIMS report
            - `depth`: the depth of coverage at the mutation site (optionally provided during initialization)
            - `frequency`: the frequency of mutation at the site (optionally provided during initialization)
            - `read_support`: the number of reads supporting the mutation (depth*frequency)
            - `rationale`: whether or not the mutation was called using the WHO confidence or expert rules
            - `warning`: a column reserved for warnings such as low depth of coverage or poor read support
        - If the originating `Variant` _is_ None, the `Row` object is initalized with the following attributes:
            - `tbprofiler_gene_name` must be provided
            - all other items (except `locus_tag`) are by default "NA", unless they pass quality metrics in which case they become "WT". If they fail these quality metrics, they become "Insufficient Coverage" and a `warning` is created (these quality metrics are measured using the `COVERAGE_DICTIONARY` created by the `Coverage` object).
    - The `Row` object (made in the `Variant.extract_annotation()` method is then appended to the `Variant`'s `annotation_dictionary` attribute which is returned to the `iterate_section` function in the `Laboratorian` object. The `iterate_section` function then iterates through each `Row` in the `annotation_dictinoary` and performs the `Row` object's `complete_row()` method.
        - The `complete_row()` method will either:
            1. If the `Row` object has a `who_confidence` attribute that is _not_ absent, it will convert that confidence into a `looker_interpretation` and `mdl_interpretation` using the `ANNOTATION_TO_INTERPRETATION` global variable dictionary.
            2. If the `Row` object is missing a `who_confidence`, it will apply the `apply_expert_rules()` method to its originating `Variant` object.
    - The `apply_expert_rules()` function in the `Variant` class implements the rules outlined in the CDPH interpretation document regarding the interpretation of potential resistance mutations. Each `Variant` is considered based on its originating gene and the antimicrobial drug it is associated with. Depending on the mutation and/or position, the appropriate Looker and MDL interpretations are returned to the `Row` object.
        - Depending on whether or not the rule applied is considered an "expert" rule or not, a `"noexpert"` string is sometimes appended to the interpretation. If this is the case, the `complete_row()` method calls the `Row`'s `describe_rationale()` function which removes this suffix and modifies the `rationale` field.
    - Following `Row` completion, depending on if the `Row`'s `tbprofiler_gene_name` belongs to a certain group of genes, the `extract_alternate_consequences` function is called on the `Variant`. This performs a function similar to the `iterate_section` method but upon a `Variant`'s `alternate_consequences` attribute. Any additional `Row`s made are added to the `row_list`.
5. After `iterate_section()` completes within the `create_laboratorian_report()` method, any genes that are not already in the `row_list` are added. In addition, any warnings are also applied. Finally, the `row_list` is added to the global variable `DF_LABORATORIAN` for usage in other functions, and the row_list is converted into an output CSV file.
6. The `Parser.run()` method then creates a `LIMS` object (`LIMS.py`), which is initialized using the input JSON as well. The `create_lims_report()` method is then called on the object.
7. The `create_lims_report()` method iterates through each `{antimicrobial: gene}` combination in the global variable `ANTIMICROBIAL_CODE_TO_GENES`. It determines the maximum MDL resistance ranking for that gene-drug combo, and then converts the annotation into LIMS language with the `convert_annotation()` method. Afterwards, the `apply_lims_rules()` function is called.
    - The `apply_lims_rules()` function iterates through each gene associated with an antimicrobial drug and formats the mutations into LIMS language. Additionally, it applies special language for "rifampin" and also determines based on position, quality scores, etc. what language should be used. Similar to `apply_expert_rules()` in the `Variant` class, this method follows the guidelines described in the CDPH interpretation document.
8. After `apply_lims_rules()` method completes, the `create_lims_report()` method saves the output to a CSV file.
9. The `Parser.run()` method now creates a `Looker` object (`Looker.py`) which is initialized with the input JSON file. The `create_looker_report()` method is then called on the object.
10. The `create_looker_report` iterates through each drug in the global variable `ANTIMICROBIAL_DRUG_NAME_LIST` and determines the highest mutation ranking for the drug; that is, the most severe resistance ranking is given as output (R > U > S). If the highest ranking is not "R" and the gene has poor coverage, the ranking is replaced with "Insufficient coverage in locus". Finally, the `create_looker_report()` method saves the output to a CSV file.
11. The `Parser.run()` method then finalizes the coverage report using the `Coverage` object's `reformat_coverage()` method, which adds a warning to every gene with a deletion and saves the output to a CSV file.
12. The `Parser.run()` method then prints a message to the user indicating that the function has completed successfully. Upon completion of this function, the program ends.
