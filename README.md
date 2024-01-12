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

```markdown
docker pull us-docker.pkg.dev/general-theiagen/theiagen/tbp-parser:1.1.9
```

## Usage

```markdown
usage: python3 /tbp-parser/tbp_parser/tbp_parser.py [-h|-v] <input_json> <input_bam> [<args>]

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
          Do not include any spaces
          default="tbp_parser"
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

### Example

```markdown
python3 /tbp-parser/tbp_parser/tbp_parser.py \
    /path/to/data/tbprofiler_output.json \
    /path/to/data/tbprofiler_output.bam \
    -o "example" \
    -d 12 \
    -c 98 \
    -s "Illumina NextSeq" \
    -p "John Doe" 
```

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

## Description of Process

What follows is a step-by-step explanation of what happens in the tbp-parser tool:

1. The `tbp_parser.py` script is run with the input JSON and BAM files. This script is the entry point for the tool, and begins by parsing the various command-line arguments. Upon completion of arguments, it creates an new `Parser` object with those arguments as input and runs the `.run()` method of that class.
2. The Parser class (`Parser.py`) is the controlling module for this package. It is initialized with the command line arguments and uses these inputs to set global variables, such as `MIN_DEPTH`, `COVERAGE_THRESHOLD`, etc. Its primary method is `run()`. This method also calls the `check_dependency_exists()` function, which determines if the prequesite `samtools` is installed. 
3. The first object created in the `run()` function is an instance of the `Coverage` class (`Coverage.py`) in order to make the initial coverage report. This class is initialized with the input BAM file. It then runs the `create_coverage_report()` method, which creates the coverage report and saves it to a CSV file.
    - The `calculate_coverage()` method iterates through the genes contained in the `../data/tbdb.bed` file (collected on Sep 1, 2023 from the TBProfiler repository). It runs `samtools depth` and then divides the sum of the positions above the `MIN_DEPTH` value by the length of the gene in question to calculate the percent coverage of each gene in the H37Rv reference genome. It then saves this information as `{gene: percent coverage}` into the global variable `COVERAGE_DICTIONARY`.
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
            - `mdl_interpretation`: the MDL interpretation of resistance
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
6. The `Parser.run()` method now creates a `Looker` object (`Looker.py`) which is initialized with the input JSON file. The `create_looker_report` method is then called on the object.
7. The `create_looker_report` iterates through each drug in the global variable `ANTIMICROBIAL_DRUG_NAME_LIST` and determines the highest mutation ranking for the drug; that is, the most severe resistance ranking is given as output (R > U > S). If the highest ranking is not "R" and the gene has poor coverage, the ranking is replaced with "Insufficient coverage in locus". Following determining the maximum ranking, the `get_lineage_and_id()` function is called which determines the lineage and ID of the sample. Finally, the `create_looker_report()` method saves the output to a CSV file.
8. The `Parser.run()` method then creates a `LIMS` object (`LIMS.py`), which is initialized using the input JSON as well. The `create_lims_report()` method is then called on the object.
9. The `create_lims_report()` method iterates through each `{antimicrobial: gene}` combination in the global variable `ANTIMICROBIAL_CODE_TO_GENES`. It determines the maximum MDL resistance ranking for that gene-drug combo, and then converts the annotation into LIMS language with the `convert_annotation()` method. Afterwards, the `apply_lims_rules()` function is called.
    - The `apply_lims_rules()` function iterates through each gene associated with an antimicrobial drug and formats the mutations into LIMS language. Additionally, it applies special language for "rifampin" and also determines based on position, quality scores, etc. what language should be used. Similar to `apply_expert_rules()` in the `Variant` class, this method follows the guidelines described in the CDPH interpretation document.
10. After `apply_lims_rules()` method completes, the `create_lims_report()` method saves the output to a CSV file.
11. The `Parser.run()` method then finalizes the coverage report using the `Coverage` object's `reformat_coverage()` method, which adds a warning to every gene with a deletion and saves the output to a CSV file.
12. The `Parser.run()` method then prints a message to the user indicating that the function has completed successfully. Upon completion of this function, the program ends.
