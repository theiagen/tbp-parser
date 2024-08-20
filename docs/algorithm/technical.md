--- 
title: "Technical Code Breakdown"
---

`tbp-parser` is object-oriented, with each class representing either *an output file*, *a part of an output file*, or *a part of the input JSON file* produced by TBProfiler.

The first class that is invoked by the `tbp-parser.py` script is `Parser` which is a control class that orchestrates the creation of the different output reports. 

## Calculating percent gene coverage

Before creating any reports, `Parser` calls the `Coverage` class to calculate the percent gene coverage over a specified minimum depth (default: 10) for the coding regions of all genes included in the TBDB (the database used in TBProfiler to generate the drug resistance annotations). This requires as input the BAM and BAI files produced by TBProfiler during alignment to the H37Rv reference genome. The percent gene coverage results are then stored in a global dictionary that is accessed multiple times for QC purposes during the creation of the final reports.

## Creating the Laboratorian report

Then, `Parser` creates the Laboratorian report using the `Laboratorian` class and its associated `.create_laboratorian_report()` method.

The `Laboratorian` class uses the input JSON file to collect the necessary information. The structure of the input JSON file is a good place to start the breakdown:

```json linenums="1" title="example_input.json"
{
  "id": "sample01",
  ...
  "main_lin": "lineage1",
  "sublin": "lineage1.2.1.2.1",
  "dr_variants": [ ... ],
  "other_variants": [ ...],
  ...
}
```

In this example, we can see only the relevant top-level JSON fields that are used in `tbp-parser`.

Of interest, the `"id"` column is used to set the global `SAMPLE_NAME` variable.

The lineage information, found in `"main_lin"` and `"sublin"` are used in the LIMS and Looker reports, so we won’t go into detail about them here.

The variant information is what makes up the bulk of the Laboratorian report and can be found in the `"dr_variants"` and `"other_variants"` fields. We’ll talk more about these fields later.

There are many other fields that are omitted from this example since they are not used in `tbp-parser`, such as version information and overall sample drug resistance type (like RR-TB, etc.). These fields are found in the other TBProfiler output files in more human-readable formats.

Within the input JSON file, there are two fields that are examined the most: `"dr_variants"` and `"other_variants"`. These fields are treated the same, and have the same format, although different mutations are found in both regions. The difference between the two fields is unclear to me at this time. In the example below, only the fields used in `tbp-parser` are shown.

```json linenums="1" title="example_input.json"
{
  ...
  "dr_variants": [
    {
      "chrom": "Chromosome",
      "genome_pos": 761109,
      ...
      "depth": 130,
      "freq": 1,
      "type": "missense_variant",
      "nucleotide_change": "c.1303G>T",
      "protein_change": "p.Asp435Tyr",
      "annotation": [
        {
          "type": "who_confidence",
          "drug": "rifampicin",
          "who_confidence": "Assoc w R"
        }
      ],
      "alternate_consequences": [],
      ...
      "locus_tag": "Rv0667",
      "gene": "rpoB",
      "gene_associated_drugs": [
        "rifampicin"
      ]
    }
  ],
  "other_variants": [
    {
      "chrom": "Chromosome",
      "genome_pos": 6112,
      ...
      "depth": 105,
      "freq": 1,
      "type": "missense_variant",
      "nucleotide_change": "c.873G>C",
      "protein_change": "p.Met291Ile",
      "annotation": [
        {
          "type": "who_confidence",
          "drug": "moxifloxacin",
          "who_confidence": "Not assoc w R"
        },
        {
          "type": "who_confidence",
          "drug": "levofloxacin",
          "who_confidence": "Not assoc w R"
        }
      ],
      "alternate_consequences": [],
      ...
      "locus_tag": "Rv0005",
      "gene": "gyrB",
      "gene_associated_drugs": [
        "levofloxacin",
        "ofloxacin",
        "moxifloxacin",
        "fluoroquinolones",
        "ciprofloxacin"
      ]
    },
...
}
```

After the global `SAMPLENAME` variable is set, the `Laboratorian` class calls the `.iterate_section()` method, starting with the `"dr_variants"` field.

Since the contents of each variant section in the JSON dictionary are considered a list, we start to iterate through each list item, which consists of each section within curly brackets `{...}`. In the example to the left, I’ve only included 1 item in each list. 

Immediately, each item in the list is converted into a `Variant` class object, and every item in each list item (the `"chrom"`, `"genome_pos"`, `"locus_tag"`, etc.) is converted to a class attribute. This is because each item in the list represents a single mutation or a single variant. I’ll now refer to each variant section item as a `Variant`.

Each new `Variant` object has the `.extract_annotations()` method called. This method starts by iterating through the `"annotation"` field in the input JSON. The annotation field can contain multiple different annotations, so we look at each one individually.

Each annotation is turned into a `Row` object, which represents a row in the Laboratorian report. During the initiation of the `Row` object, each column in the Laboratorian report is created based on both the annotation field and the originating `Variant` object. Additionally, a warning field is created based on both the global dictionary created with the `Coverage` class and the mutation’s `"depth"` and `"freq"` fields.

Sometimes multiple annotations for the same drug can appear for a single `Variant`. If this is the case, only the most severe annotation is saved (that is, an annotation that indicates resistance is kept instead of one that indicates susceptibility).

After the annotation field has been iterated through, we then check the `"gene_associated_drugs"` field to make sure that we create a `Row` for each antimicrobial drug that is associated with the gene. As you can see in the `"other_variants"` section, the annotation field for the variant only lists annotations for moxifloxacin and levofloxacin, but the gene is associated with three other antimicrobial drugs. This iteration creates additional `Row` objects for those antimicrobial drugs.

This means that each mutation will potentially appear several times in the final report, once for every antimicrobial associated with the drug. This is because sometimes a mutation confers a different resistance level to one drug, but not another.

After `Row` objects are created for each `Variant` in the variant section, every `Row` has the `.complete_row()` method called, which adds the interpretation columns to the object. Two interpretation columns are created, `mdl_interpretation` and `looker_interpretation`.

Please note that these interpretation columns are typically identical, but in several cases, the `mdl_interpretation` column will call a variant-drug combination as “susceptible” (S), while the `looker_interpretation` column will call the same combination “uncertain” (U).

In the case where a WHO annotation was not identified, the `Variant` class’ `.apply_expert_rules()` method is called. This function applies expert rules that are listed in detail on the `tbp-parser` GitHub repository, [available here](https://github.com/theiagen/tbp-parser/blob/49407cea8088886d76d772353f647776f9280d7c/2024-01-30-interpretation_logic.pdf).

The expert rules assign a drug resistance call to the variant-drug combination only when there is no WHO annotation and will fill the `mdl_interpretation` and `looker_interpretation` fields.

If the mutation is in either mmpS5, mmpL5, or mmpR5/Rv0678, then the `"alternate_consequences"` field is iterated through. This field typically lists the same mutation but in reference to a different gene; for instance, if a mutation is in the upstream non-coding region of one gene, it may be in the coding region of a different gene.

Then, any genes that do not have any variants are added to the laboratorian report with various “NA” or “WT” values filling the appropriate fields.

This means that every gene in the TBDB appears in the Laboratorian report regardless if any mutations were identified in that gene.

Finally, a few more quality control measures are taken and then all of the individual `Row` objects are written to a CSV file, which concludes the creation of the laboratorian report.

## Creating the Looker report

The `Parser` class then creates a `Looker` object which uses the `.create_looker_report()` method. The Looker report uses the Laboratorian report to generate most of the included information.

It starts by iterating through a list of antimicrobial drugs and extracting all of the `looker_interpretation` values for each row in the report with that antimicrobial drug. It then identifies the highest resistance rating (R > R-Interim > U > S-Interim > S) for all resistance annotations for a drug.

Then, a quality check is performed and if a particular gene fails coverage that contributed to the highest resistance rating, an insufficient coverage warning is given.

The `"main_lin"` and `"sublin"` fields from the input JSON file are used to fill the `ID` field in the report. These fields are converted into shortened English without any technical lineage information.

Finally, the information is written to a CSV file which concludes the creation of the Looker report.

## Creating the LIMS report

The `Parser` class then creates `LIMS` object which uses the `.create_lims_report()` method. The LIMS report also uses the Laboratorian report to generate the bulk of the information included.

The `.create_lims_report()` method begins by iterating through each LIMS antimicrobial and gene code (corresponding to the LIMS codes in the CDPH STAR LIMS system). Then, the highest `mdl_interpretation` value is extracted for each row in the report that is associated with that antimicrobial drug, like in the Looker report. Then, the annotation is converted into a human-readable format (R → Mutations(s) associated with resistance to {antimicrobial} detected”, etc.).

Then, the `.apply_lims_rules()` function is activated which determines which mutations should be output for the corresponding drug-gene combination. The mutations are then formatted so that they appear in the following format: `{nucleotide mutation} ({amino acid mutation, if available})` repeated, separated by semicolons.

Some specific parsing rules apply to mutations within the rpoB gene, which changes the output language on the LIMS report. These rules depend on the position of the mutation in the gene.

After the rules are applied and the mutations are collected, the information is written to a CSV file which concludes the creation of the LIMS report.

## Creating the coverage report

The `Parser` class then reuses the `Coverage` object created first and calls the `.reformat_coverage()` method which adds any warnings, such as any deletion mutations detected for a gene. If a deletion is detected, a warning is useful because it indicates that although the reported coverage is less than 100%, it may be due to that deletion. If the coverage is still 100% and a deletion was identified, the warning will say that the deletion *may* be upstream.

The coverage dictionary and the associated warnings are then written to a CSV file which concludes the creation of the coverage report, and the `tbp-parser` script.
