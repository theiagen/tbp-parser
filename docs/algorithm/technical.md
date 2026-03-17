# Technical Code Breakdown

With the v3.0.0 refactor, `tbp-parser` has undergone significant changes and the previous technical description no longer applies. This document serves as a technical breakdown of the new code structure and algorithm used in `tbp-parser` v3.0.0.

## Class descriptions

_In order of appearance_

- `Configuration`: a singleton class that handles all configuration settings and input parameters for `tbp-parser`
- `GeneDatabase`: a singleton class that contains information regarding each gene of interest; attributes include gene name, locus tag, tier, associated drugs, and promoter regions.
- `LIMSGeneCode`: a class representing a gene-specific result for the LIMS report 
- `LIMSRecord`: a class representing a drug-specific result for the LIMS report; contains a list of `LIMSGeneCode` objects (the genes associated with that drug)
- `BedRecord`: a class representing a record/entry from a BED file; contains attributes such as chromosome, start, end, and gene name, It derives length and genomic coordinates.
- `CoverageCalculator`: a class used to generate coverage statistics for BED records; contains methods to calculate percent coverage and average depth for a given BED record based on the read depth information from the BAM file
- `BaseCoverage`: a base class to represent coverage data fora  genomic region that can be extended to represent coverage for specific types of regions.
- `TargetCoverage`: a class to represent coverage data for a single target genomic region (BedRecord); based on the unique gene_name in the BedRecord. Multiple TargetCoveage objects may exist for a single locus_tag/gene if there are multiple BedRecords for the same locus_tag/gene.
- `LocusCoverage`: a class to represent coverage data for a single locus tag, aggregating multiple target regions if necessary.
- `ERRCoverage`: a class to represent coverage data for the essential reportable range (ERR) of a single target/locus
- `VariantRecord`: a data class represnting a single variant record (aka JSON entry) from the TBProfiler JSON output
- `VariantProcessor`: a class that processes VariantRecords into Variant objects
- `Variant`: a class that represents a single genetic drug-gene-variant combination
- `Helper`: a class containing a variety of helper functions for parsing and analyzing mutations
- `Consequences`: a data class representing a single consequence entry (dictionary) for a list of dicts under the `consequences` field in the TBProfiler JSON output
- `Annotation`: a data class representing a single annotation entry (dictionary) for a list of dicts under the `annotation` field in the TBProfiler JSON output
- `VariantInterpreter`: a class to handle interpretation of variants and expert rules
- `InterpretationResult`: a class to represent the result of variant interpretation
- `VariantQC`: a class to generate QC warnings related to variants, such as low coverage or poor quality
- `QCResult`: a data class representing the result of a QC check on a variant
- `LIMSProcessor`: a class to handle the decision logic for LIMS report generation

## Walkthrough

### Setup and input parsing

Arguments are parsed into the singleton instance of the `Configuration` class, which is used to store all input parameters. Alternatively, if a [configuration file](../inputs.md#configuration-file) is provided, it will be used to overwrite any provided command-line arguments. 

The provided (or default) [gene database YAML file](../inputs.md#gene-database-file) is then parsed into the singleton instance of the `GeneDatabase` class, which contains all relevant information for each gene of interest. This class takes the form of a dictionary, where the keys are the locus tags of the genes and the values are dictionaries containing all relevant information for that gene: gene name, locus tag, tier, associated drugs, and promoter regions. This database also exists in an "inverted" form where the key is the gene name instead of the locus tag.

The provided (or default) [LIMS report format YAML file](../inputs.md#lims-report-format-yaml-file) is then parsed into a `LIMSRecord` list and any corresponding `LIMSGeneCode` list. See the toggle below for more details.

??? techdetails "The `LIMSRecord` and `LIMSGeneCode` classes"
    `LIMSRecord` is a drug-specific result for the LIMS report. It contains the following attributes:
        
    - `drug`: the name of the drug (as it would appear in the input TBProfiler JSON file under the `"annotation.drug"` field; e.g., "bedaquiline")
    - `drug_code`: the desired output column name for that drug in the LIMS report (e.g., "BDQ")
    - `gene_codes`: a dictionary of gene names (as they appear in TBProfiler under the `gene_name` field; e.g., "mmpR5") to a corresponding `LIMSGeneCode` object
    
    `LIMSGeneCode` is a result for a gene-drug combination for the LIMS report with the following attributes:
        
      - `gene_code`: the desired output column name for that gene-drug combination in the LIMS report (e.g., "BDQ_Rv0678")
      - `gene_target_value`: the value that will be reported in the LIMS report for that gene-drug combination; see the `LIMSProcessor` class for more details on how this is determined
      - `max_mdl_interpretations`: a list of the highest `mdl_interpretation` values identified for any mutation in that gene for that drug; this is used in the decision logic for determining the `gene_target_value`
      - `max_mdl_variants`: a list of the specific mutations that have the highest `mdl_interpretation` values for that gene-drug combination; this is used in the decision logic for determining the `gene_target_value`

The provided (or default) [coverage BED file](../inputs.md#coverage-bed-file) is then parsed into a `BedRecord` list, where each `BedRecord` represents a single record/entry from the BED file., See the toggle below for more details. This same process is applied to the [ERR coverage BED file](../inputs.md#tngs-specific-arguments) (if provided), resulting in a `BedRecord` list representing the ERR regions.

??? techdetails "The `BedRecord` class"
    `BedRecord` is a class representing a single record/entry from a BED file. It contains the following attributes:

    - `chrom`: the chromosome of the region (e.g., "chromosome")
    - `start`: the start position of the region (e.g., 100)
    - `end`: the end position of the region (e.g., 200)
    - `locus_tag`: the locus tag of the gene associated with that region (e.g., "Rv0678")
    - `gene_name`: the gene name of the gene associated with that region (e.g., "mmpR5")
    - `length`: the length of the region, derived from the start and end positions (e.g., 100)
    - `coords`:  a tuple of the coordinates of the region, derived from the start and end positions (e.g., (100, 200))
    - `reads_by_position`: a dictionary where the keys are genomic positions and the values are lists of read names that appear at that position; this is populated during the coverage calculation step (see below)

Once the input database files have been processed and parsed, the identified `LIMSRecord` and `BedRecord` lists are compared to ensure that each `LIMSRecord` gene has a corresponding `BedRecord`. If any genes are missing from the BED file, an error is raised. Additionally, an error is raised if any genes in the `LIMSRecord` list are not found in the `GeneDatabase`. This ensures that all genes of interest have all information necessary for downstream processing and report generation.

### Coverage calculations

/// html | div[style='float: left; width: 50%; padding: 20px;']

After initial setup and input parameter processing, the `CoverageCalculator` class is initialized using the `Configuration` instance, which provides access to the BAM file and other necessary parameters for the breadth of coverage calculation. The `.calculate()` method is called on the `CoverageCalculator` object, which iterates through the `BedRecord` list (and the ERR `BedRecord` list if applicable). 

Each record in the list is run through `pysam AlignmentFile.pileup()` to identify which reads in the BAM file are associated with that record's region. Those reads are stored in the `BedRecord`'s `"reads_by_position"` attribute. Each `BedRecord` has a unique `reads_by_position` dictionary, so if regions overlap, the same positions and reads may be associated with multiple `BedRecord` entries.

///

/// html | div[style='float: right; width: 50%; padding: 20px;']

!!! techdetails "`reads_by_position`"
    `reads_by_position` is a dictionary that takes the following format, where the keys are genomic positions and the values are lists of read names that appear at that position:

    ```python
    {
        100: [read1, read2, read3],
        101: [read2, read3],
        102: [read3, read4],
        ...
    }
    ```

///

/// html | div[style='clear: both;']
///

??? techdetails "Handling overlapping primer regions in tNGS analyses"
    If `--tNGS` is used, the `BedRecord` list is checked to determine if any two `BedRecord` coordinates overlap. This is done by whitelisting any non-overlapping reads for each `BedRecord` under the assumption that if a read appears in a non-overlapping region it originated from that particular target. The `reads_by_position` attribute is then filtered to only include those whitelisted reads. 
    
    For example, imagine geneA is covered by two primers: primer1 covers bases 0-45, and primer2 covers bases 30-75, meaning that 15 bases overlap between the primers. readA appears in the `reads_by_position` dictionary for primerA positions 0-45, and appears in the primerB `reads_by_position` dictionary for positions 30-45. readB appears in the primerB dictionary for positions 30-75 and the primerA dictionary for positions 30-45. readA appears in the **non-overlapping** region of primer1 (0-30) while readB appears in the **non-overlapping region** of primer2 (45-75). readA is then **whitelisted** for primer1 and readB is **whitelisted** for primer2. 
    
    If a third read, readC, covers only positions 30-45, this means that it appears only in the overlapping region of both primer1 and primer2. It is not whitelisted for either record, and so it is filtered out of the `reads_by_position` attribute for both records since its origin cannot be determined. 
    
    This process allows us to make an informed assumption about which reads originated from which target regions, which is important for accurate coverage calculations of individual target regions in the tNGS analysis, since it isolates the reads belonging to each primer.

Once the `BedRecord` list is finalized, the `reads_by_position` dictionary is used to calculate breadth of coverage. The number of positions that has more reads than the number specified by `--min_depth` is divided by the number of positions in the dictionary. Average depth is calculated by summing the number of reads at each position and then dividing by the length of the dictionary. These coverage statistics are stored in the `BedRecord`'s `"breadth_of_coverage"` and `average_depth` attributes.

If an ERR coverage BED file is provided, the same process is applied to the ERR `BedRecord` list, resulting in breadth of coverage and average depth calculations for the ERR regions. These regions are then associated with their corresponding BRR regions based on the gene name and locus tag under `err_coverage`.

### TBProfiler JSON processing

The input TBProfiler JSON is then parsed. Included are examples of the input JSON format, with explanations of the relevant fields for `tbp-parser` and how they are used in the algorithm.

/// html | div[style='float: left; width: 50%; padding: 20px;']

```json linenums="1" title="example_input.json: top level fields"
{
    "schema_version": "1.0.0",
    "id": "sample01",
    "timestamp": "2025-04-25T22:12:11.677658",
    "pipeline": { ... },
    "notes": [],
    "lineage": [ ... ],
    "main_lineage": "lineage4",
    "sub_lineage": "lineage4.3.4.1",
    "spoligotype": null,
    "drtype": "XDR-TB",
    "dr_variants": [ ... ],
    "other_variants": [ ...],
    "qc_fail_variants": [],
    "qc": { ... },
    "linked_samples": []
}
```

///

/// html | div[style='float: right; width: 50%; padding: 20px;']

In the example to the left, we see only the top-level JSON fields, only some of which are used in tbp-parser.

The `"id"` field is saved to the `SAMPLE_ID` variable, and the lineage information, found in the `"main_lineage"` and `"sub_lineage"` fields, is extracted and saved in the `LINEAGE_ID` and `SUBLINEAGE_ID` variables respectively.

The variant information is what makes up the bulk of the Laboratorian report and can be found in the `"dr_variants"` and `"other_variants"` fields. Each of these fields contains a list of variants, where each variant is represented by a dictionary with many different fields. The relevant fields for `tbp-parser` are explained in the next section.

There are many other fields available that have useful information but are not used in tbp-parser, such as versioning (found in `"pipeline"`), quality control metrics (found in `"qc"`) and overall sample drug resistance type (found in `"drtype"`, like RR-TB, etc.). These fields can also be found in the other TBProfiler output files in more human-readable formats.

///

/// html | div[style='clear: both;']
///

/// html | div[style='float: right; width: 50%; padding: 20px;']

```json linenums="1" title="example_input.json: dr_variants and other_variants"
{
  ...
  "dr_variants": [
    {
          "chrom": "Chromosome",
          "pos": 779127,
          "ref": "T",
          "alt": "TG",
          "depth": 109,
          "freq": 0.5137614678899083,
          "sv": false,
          "filter": "pass",
          "forward_reads": 24,
          "reverse_reads": 32,
          "sv_len": null,
          "gene_id": "Rv0678",
          "gene_name": "mmpR5",
          "feature_id": "CCP43421",
          "type": "frameshift_variant",
          "change": "c.139dupG",
          "nucleotide_change": "c.139dupG",
          "protein_change": "p.Asp47fs",
          "annotation": [ ... ],
          "consequences": [ ... ],
          "drugs": [ ... ],
          "locus_tag": "Rv0678",
          "gene_associated_drugs": [
              "bedaquiline",
              "clofazimine"
          ]
      },
      ...
  ],
  "other_variants": [
      { 
          "chrom": "Chromosome",
          "pos": 3065027,
          "ref": "G",
          "alt": "A",
          "depth": 94,
          "freq": 1.0,
          "sv": false,
          "filter": "pass",
          "forward_reads": 41,
          "reverse_reads": 53,
          "sv_len": null,
          "gene_id": "Rv2752c",
          "gene_name": "Rv2752c",
          "feature_id": "CCP45551",
          "type": "missense_variant",
          "change": "p.Arg389Trp",
          "nucleotide_change": "c.1165C>T",
          "protein_change": "p.Arg389Trp",
          "annotation": [ ... ],
          "consequences": [ ... ],
          "locus_tag": "Rv2752c",
          "gene_associated_drugs": [
                "ethambutol",
                "isoniazid",
                "levofloxacin",
                "moxifloxacin",
                "rifampicin"
          ]
      },
      ...
  ]
...
}
```

///

/// html | div[style='float: left; width: 50%; padding: 20px;']

Each entry in `"dr_variants"` and `"other_variants"` represents a single variant identified by TBProfiler. Each item is saved to a `VariantRecord` object, which is a data class representing a single variant record (aka JSON entry) from the TBProfiler JSON output, containing the following attributes:

!!! techdetails "`VariantRecord` attributes"
    - `sample_id`: the sample ID (top-level `"id"`)
    - `pos`: the genomic position of the variant (`pos`)
    - `depth`: the depth of coverage at the variant position (`depth`)
    - `freq`: the frequency of the variant in the reads (`freq`)
    - `gene_id`: the locus tag of the gene associated with the variant (`gene_id`)
    - `gene_name`: the gene name associated with the variant (`gene_name`)
    - `type`: the type of variant (e.g., frameshift, missense, etc.) (`type`)
    - `nucleotide_change`: the nucleotide change associated with the variant (`nucleotide_change`)
    - `protein_change`: the protein change associated with the variant (`protein_change`)
    - `annotation`: a list of `Annotation` objects associated with the variant (`annotation`)
    - `consequences`: a list of `Consequence` objects associated with the variant (`consequences`)
    - `gene_associated_drugs`: a list of drugs associated with the gene of the variant (`gene_associated_drugs`)

///

/// html | div[style='clear: both;']
///

All `VariantRecord` objects are added to a list which is provided to the `VariantProcessor` class for downstream processing, including expansion, annotation extraction, interpretation, and QC. This method class begins by processing the `VariantRecord` list by first, expanding the `"consequences"` attribute of the `VariantRecord` into additional `VariantRecord` objects. This occurs only if the `VariantRecord` has a `gene_id` of Rv0676c (_mmpL5_), Rv0677c (_mmpS5_), or Rv0678 (_mmpR5_). 

/// html | div[style='float: left; width: 50%; padding: 20px;']

Each entry in the `"consequences"` list is a dictionary, and each of those dictionaries is saved as an attribute of a `Consequence` object. If this field is blank or missing, an empty `Consequence` object is still created to enable ease of downstream processing.

!!! techdetails "`Consequence` attributes"
    - `gene_id`: the locus tag of the gene associated with the consequence (e.g., "Rv0676c")
    - `gene_name`: the gene name associated with the consequence (e.g., "mmpL5")
    - `type`: the type of consequence (e.g., "upstream_gene_variant")
    - `nucleotide_change`: the nucleotide change associated with the consequence (e.g., "c.-648dupC")
    - `protein_change`: the protein change associated with the consequence, if applicable (e.g., "")
    - `annotation`: a list of `Annotation` objects associated with the consequence (`annotation`)

Each consequence represents an alternate mapping of the same variant to a different gene. In the example above, the same variant is mapped to both mmpR5 and mmpL5, with different consequence types and annotations for each gene. This allows us to capture all possible drug resistance implications of a given variant, even if it maps to multiple genes.

///

/// html | div[style='float: right; width: 50%; padding: 20px;']

```json linenums="1" title="example_input.json: consequences field"
"consequences": [
    {
        "gene_id": "Rv0678",
        "gene_name": "mmpR5", // this entry is ignored because this is the same gene as the parent VariantRecord
        "feature_id": "CCP43421",
        "type": "frameshift_variant",
        "nucleotide_change": "c.139dupG",
        "protein_change": "p.Asp47fs",
        "annotation": [ ... ], 
    },
    {
        "gene_id": "Rv0676c",
        "gene_name": "mmpL5",
        "feature_id": "CCP43419",
        "type": "upstream_gene_variant",
        "nucleotide_change": "c.-648dupC",
        "protein_change": "",
        "annotation": []
    },
    {
        "gene_id": "Rv0677c",
        "gene_name": "mmpS5",
        "feature_id": "CCP43420",
        "type": "upstream_gene_variant",
        "nucleotide_change": "c.-223dupC",
        "protein_change": "",
        "annotation": []
    }
]
```

///

/// html | div[style='clear: both;']
///


The `_expand_consequences` method will create a copy of the variant record for each entry in the `consequences` list, replacing the gene information (gene_id, gene_name, etc.) with the corresponding information from the consequence entry. The consequence annotations are also added to the `annotation` field of the new `VariantRecord`. If the `annotation` field is blank, `Annotation` object(s) will be created that match the parent `VariantRecord` drug(s) with a `"confidence"` of "No WHO annotation" and blank `"source"` and `"comment"` fields.

/// html | div[style='float: left; width: 50%; padding: 20px;']

```json linenums="1" title="example_input.json: annotation field"
"annotation": [
    {
        "type": "drug_resistance",
        "drug": "bedaquiline",
        "original_mutation": "LoF",
        "confidence": "Assoc w R",
        "source": "WHO catalogue v2",
        "comment": "Can only confer resistance if genetically linked to a functional MmpL5"
    },
    {
        "type": "drug_resistance",
        "drug": "clofazimine",
        "original_mutation": "LoF",
        "confidence": "Assoc w R",
        "source": "WHO catalogue v2",
        "comment": "Can only confer resistance if genetically linked to a functional MmpL5"
    }
],
```

///

/// html | div[style='float: right; width: 50%; padding: 20px;']

In the example to the right, we see the `"annotation"` field for a single variant. This field contains a list of annotations, where each annotation is represented by a dictionary with many different fields. Each field in the annotation dictionary is saved as an attribute of an `Annotation` object (if the field is present and has content).

!!! techdetails "`Annotation` attributes"
    - `drug`: the drug associated with the annotation (e.g., "bedaquiline")
    - `confidence`: the WHO confidence level of the annotation (e.g., "Assoc w R")
    - `source`: the source of the annotation (e.g., "WHO catalogue v2")
    - `comment`: any comments associated with the annotation (e.g., "Can only confer resistance if genetically linked to a functional MmpL5")

In this example, this variant confers resistance to both bedaquiline and clofazimine. Each annotation is saved to the `VariantRecord` to enable ease of downstream processing. 

An annotation can be found in either a VariantRecord or in a Consequence object, since both variants and consequences can have unique annotations.

///

/// html | div[style='clear: both;']
///


