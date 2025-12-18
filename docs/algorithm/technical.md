# Technical Code Breakdown

`tbp-parser` is object-oriented, with each class representing either _an output file_, _a part of an output file_, or _a part of the input JSON file_ produced by TBProfiler.

The first class that is invoked by the `tbp-parser.py` script is `Parser` which is a control class that orchestrates the creation of the different output reports.

## Calculating percent gene coverage

Before creating any reports, `Parser` calls the `Coverage` class to calculate the percent breadth of coverage of a particular gene or locus over a specified minimum depth (default: 10) for any regions included in the `--tbdb_bed` input file (default value is the `genes.bed` file used in TBProfiler's TBDB). This requires the BAM and BAI files produced by TBProfiler during alignment to the H37Rv reference genome. The percent breadth of coverage results are then stored in a dictionary that is accessed multiple times for QC purposes during the creation of the final reports. At the same time, the average locus depth is calculated and any loci that do not meet the minimum breadth of coverage (determined by `--min_percent_coverage`, default: 1.0) are added to a list of failed loci for later use in QC checks.

## Creating the Laboratorian report

Then, `Parser` creates the Laboratorian report using the `Laboratorian` class and its associated `.create_laboratorian_report()` method.

The `Laboratorian` class uses the input JSON file to collect the necessary information. The structure of the input JSON file is a good place to start the breakdown:

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

In the example to the left, we see only the top-level JSON fields, only some of which are used in `tbp-parser`.

Of interest, the `"id"` column is used to set the global `SAMPLE_NAME` variable.

The lineage information, found in `"main_lineage"` and `"sub_lineage"` are used in the LIMS and Looker reports, so we won’t go into detail about them here.

The variant information is what makes up the bulk of the Laboratorian report and can be found in the `"dr_variants"` and `"other_variants"` fields. We’ll talk more about these fields later.

There are many other fields available that have useful information but are not used in `tbp-parser`, such as versioning (found in `"pipeline"`), quality control metrics (found in `"qc"`) and overall sample drug resistance type (found in `"drtype"`, like RR-TB, etc.). These fields can also be found in the other TBProfiler output files in more human-readable formats.

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

Within the input JSON file, there are two fields that are examined the most by tbp-parser: `"dr_variants"` and `"other_variants"`. These fields are treated the same, and have the same format, although different mutations are found in both regions.

After the sample name is retrieved from the top-level `"id"` field, the `Laboratorian` class calls the `.iterate_section()` method, which starts by iterating through the `"dr_variants"` field.

The contents of each variant section in the JSON dictionary are considered a list. Each list item consists of the contents within the curly brackets `{...}`. In the example to the right, I’ve only included 1 item in the `"dr_variants"` and `"other_variants"` list for simplicity, with all of the internal list items collapsed for simplicity.

Each item in each of the two lists is converted into a `Variant` class object, and every item in each list item (the `"chrom"`, `"pos"`, `"locus_tag"`, etc.) is converted to a class attribute. I’ll now refer to each variant section item as a `Variant`.

Each new `Variant` object has the `.extract_annotations()` method called. This method starts by iterating through the `"annotation"` field in the input JSON. The annotation field can contain multiple different annotations, so we look at each one individually.

///

/// html | div[style='clear: both;']
///

/// html | div[style='float: left; width: 50%; padding: 20px;']

```json linenums="1" title="example_input.json: annotation field"
{
    ...
    "dr_variants": [
        {
            ...
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
            ...
        },
        ...
    ]
    ...
}
```

///

/// html | div[style='float: right; width: 50%; padding: 20px;']

Each item in the `"annotation"` is turned into a `Row` object, which represents a row in the Laboratorian report. During the initiation of the `Row` object, each column in the Laboratorian report is created based on both the annotation field and the originating `Variant` object. Warnings are determined based on the results from the `Coverage` calculations and the `Variant` object's `"depth"` and `"freq"` fields.

Sometimes multiple annotations for the same drug can appear for a single `Variant`. If this is the case, only the most severe annotation is saved (that is, an annotation that indicates resistance is kept instead of one that indicates susceptibility).

After the annotation field has been iterated through, we then check the `"gene_associated_drugs"` field to make sure that we create a `Row` for each antimicrobial drug that is associated with the gene. For example, if the annotation field is missing some of the drugs that are present in the `"gene_associated_drugs"` field, this iteration creates additional `Row` objects for those missing antimicrobial drugs.

This means that each mutation will potentially appear several times in the final report, once for every antimicrobial associated with the drug. This is because sometimes a mutation confers a different resistance level to one drug, but not another.

///

/// html | div[style='clear: both;']
///

After `Row` objects are created for each `Variant` in the variant section, every `Row` has the `.complete_row()` method called, which adds the interpretation columns to the object. Two interpretation columns are created, `mdl_interpretation` and `looker_interpretation`.

Please note that these interpretation columns are typically identical, but in several cases, the `mdl_interpretation` column will call a variant-drug combination as “susceptible” (S), while the `looker_interpretation` column will call the same combination “uncertain” (U).

In the case where a WHO annotation was not identified, the `Variant` class’ `.apply_expert_rules()` method is called. This function applies expert rules that are listed in detail on the `tbp-parser` GitHub repository in the interpretation document.

The expert rules assign a drug resistance call to the variant-drug combination only when there is no WHO annotation and will fill the `mdl_interpretation` and `looker_interpretation` fields.

/// html | div[style='float: right; width: 50%; padding: 20px;']

```json linenums="1" title="example_input.json: consequences field"
{
    ...
    "dr_variants": [
        {
            ...
            "consequences": [
                {
                    "gene_id": "Rv0678",
                    "gene_name": "mmpR5",
                    "feature_id": "CCP43421",
                    "type": "frameshift_variant",
                    "nucleotide_change": "c.139dupG",
                    "protein_change": "p.Asp47fs",
                    "annotation": [ ... ], // ignored because this is the same gene as the originating Variant
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
            ],
            ...
        },
        ...
    ]
    ...
}
```

///

/// html | div[style='float: left; width: 50%; padding: 20px;']

If the mutation is in either mmpS5, mmpL5, or mmpR5/Rv0678, then the `"consequences"` field is iterated through. This field typically lists the same mutation with relation to the same gene, but can occasionally list that mutation in reference to a different gene; for instance, if a mutation is in the upstream non-coding region of one gene, it may be in the coding region of a different gene.

Each consequence that is _not_ in the same gene as the originating `Variant` is added to a new `Row` object. The only items inherited from the originating `Variant` are the `"chrom"`, `"pos"`, `"depth"`, and `"freq"` fields. The rest of the fields are filled in with the information from the consequence.

The annotation fields of each consequence are treated as original annotations, and the same process is followed as described above. This means that each consequence can have its own set of annotations, which may or may not be the same as the originating `Variant`.

///

/// html | div[style='clear: both;']
///

Any genes that do not have any variants are added to the laboratorian report with various "NA" or "WT" values filling the appropriate fields.

This means that every gene in the TBDB appears in the Laboratorian report regardless if any mutations were identified in that gene.

Finally, a few more quality control measures are taken and then all of the individual `Row` objects are written to a CSV file, which concludes the creation of the laboratorian report.

Here's a breakdown of where every column in the Laboratorian report that is directly extracted from the TBProfiler output JSON. Please note that when a `Variant` object is mentioned, it corresponds to either an entry in the `"dr_variants"` or `"other_variants"` field, or a `"consequence"` of one of those entries.

| Column Name | Source in JSON file |
| --- | --- |
| sample_id | the top-level `id` field |
| tbprofiler_gene_name | the `gene_name` field in a `Variant` object |
| tbprofiler_locus_tag | the `locus_tag` field in a `Variant` object |
| tbprofiler_variant_substitution_type | the `type` field in a `Variant` object |
| tbprofiler_variant_substitution_nt | the `nucleotide_change` field in a `Variant` object |
| tbprofiler_variant_substitution_aa | the `protein_change` field in a `Variant` object |
| confidence | the `confidence` field from the `annotation` field in a `Variant` object (or manually written "No WHO annotation" if `annotation` has no content) |
| antimicrobial | the `drug` field from the `annotation` field in a `Variant` object (or a drug in the `gene_associated_drugs` list if `annotation` has no content) |
| depth | the `depth` field in a `Variant` object |
| frequency | the `freq` field in a `Variant` object |
| source | the `source` field from the `annotation` field in a `Variant` object (or left blank if `annotation` has no content) |
| comment | the `comment` field from the `annotation` field in a `Variant` object (or left blank if `annotation` has no content) |

## Creating the Looker report

The `Parser` class then creates a `Looker` object which uses the `.create_looker_report()` method. The Looker report uses the Laboratorian report to generate most of the included information.

It starts by iterating through a list of antimicrobial drugs and extracting all of the `looker_interpretation` values for each row in the report with that antimicrobial drug. It then identifies the highest resistance rating (R > R-Interim > U > S-Interim > S) for all resistance annotations for a drug.

Then, a quality check is performed and if a particular gene fails coverage that contributed to the highest resistance rating, an insufficient coverage warning is given.

The `"main_lineage"` and `"sub_lineage"` fields from the input JSON file are used to fill the `ID` field in the report. These fields are converted into shortened English without any technical lineage information.

Finally, the information is written to a CSV file which concludes the creation of the Looker report.

## Creating the LIMS report

The `Parser` class then creates `LIMS` object which uses the `.create_lims_report()` method. The LIMS report also uses the Laboratorian report to generate the bulk of the information included.

The `.create_lims_report()` method begins by iterating through each LIMS antimicrobial and gene code (corresponding to the LIMS codes in the CDPH STAR LIMS system). Then, the highest `mdl_interpretation` value is extracted for each row in the report that is associated with that antimicrobial drug, like in the Looker report. Then, the annotation is converted into a human-readable format (R → "Mutations(s) associated with resistance to {antimicrobial} detected", etc.).

Some specific parsing rules apply to mutations within the rpoB gene, which changes the output language on the LIMS report. These rules depend on the position of the mutation in the gene.

After the rules are applied and the mutations are collected, the information is written to a CSV file which concludes the creation of the LIMS report.

## Creating the coverage report

The `Parser` class then reuses the `Coverage` object created first and calls the `.create_coverage_report()` method which adds any warnings, such as any deletion mutations detected for a gene. If a deletion is detected, a warning is useful because it indicates that although the reported coverage is less than 100%, it may be due to that deletion. If the coverage is still 100% and a deletion was identified, the warning will say that the deletion _may_ be upstream.

The coverage dictionary and the associated warnings are then written to a CSV file which concludes the creation of the coverage report, and the `tbp-parser` script.
