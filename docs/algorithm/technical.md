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


