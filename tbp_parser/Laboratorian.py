from Row import Row
from Variant import Variant
import globals as globals_
import pandas as pd
import json

class Laboratorian:
    """
    This class creates the CDPH Laboratorian report.
    It has two functions:
        - iterate_section: iterates through each variant in a section and creates
            a Variant in the Variant class for each variant
        - create_laboratorian_report: creates the laboratorian report CSV file
    """

    def __init__(self, logger, input_json, output_prefix, min_depth, min_frequency, min_read_support, coverage_dictionary, low_depth_of_coverage_list, gene_to_antimicrobial_drug_name, gene_to_locus_tag, tngs_regions, gene_to_tier, promoter_regions, tngs):
        self.logger = logger
        self.input_json = input_json
        self.output_prefix = output_prefix
        
        # quality control variables
        self.MIN_DEPTH = min_depth
        """An integer representing the minimum depth of coverage required for QC."""
        self.MIN_FREQUENCY = min_frequency
        """A float representing the minimum frequency required for QC."""
        self.MIN_READ_SUPPORT = min_read_support
        """An integer representing the minimum read support required for QC."""
        
        # look-up dictionaries and lists 
        self.COVERAGE_DICTIONARY = coverage_dictionary
        """A dictionary containing coverage information for each gene/locus."""
        self.LOW_DEPTH_OF_COVERAGE_LIST = low_depth_of_coverage_list
        """A list of genes/loci that have low depth of coverage."""
        self.GENE_TO_ANTIMICROBIAL_DRUG_NAME = gene_to_antimicrobial_drug_name
        """A dictionary mapping genes to their associated antimicrobial drug names."""
        self.GENE_TO_LOCUS_TAG = gene_to_locus_tag
        """A dictionary mapping genes to their locus tags."""
        self.TNGS_REGIONS = tngs_regions
        """A dictionary containing the specific tNGS regions for each gene (to account for non-overlapping primers in the same gene)."""
        self.GENE_TO_TIER = gene_to_tier
        """A dictionary mapping genes to their tiers."""
        self.PROMOTER_REGIONS = promoter_regions
        """A dictionary containing promoter regions for each gene."""
        self.TNGS = tngs
        """A boolean indicating whether tNGS mode is enabled."""

        self.genes_reported = set()
        """A set of genes that have already been reported in the laboratorian report."""
        self.genes_with_valid_deletions = set()
        """A set of genes that have valid deletions reported in the laboratorian report."""
        self.positional_qc_fails = set()
        """A set of mutations that have failed POSITIONAL QC checks."""

    def iterate_section(self, variant_section, row_list) -> list[Row]:
        """This function iterates through each subsection in the "dr_variants" and "other_variants" sections
        in the TBProfiler JSON file. Each item in the subsection is converted into an individual Variant object, 
        and each annotation within that Variant is extracted and converted into a Row object, where each Row is
        a row in the Laboratorian report.

        Args:
            variant_section (json.load() object): The subsection of the TBProfiler JSON file to iterate through
            row_list (List): A list of rows that have been extracted so far

        Returns:
            list[Row]: The list of rows, including the newly extracted rows
        """        
        self.logger.debug("LAB:iterate_section:Iterating through the variant section to turn each one into a Variant object")
        for variant in variant_section:
            # create a Variant object and add the origin gene to the class genes_reported set variable
            variant = Variant(self.logger, self.SAMPLE_NAME, 
                              self.GENE_TO_LOCUS_TAG, self.GENE_TO_TIER, self.GENE_TO_ANTIMICROBIAL_DRUG_NAME, self.PROMOTER_REGIONS, variant)
            self.genes_reported.add(variant.gene_name)

            # tNGS only: renaming the gene to the segment name to get the coverage for QC
            if self.TNGS and variant.gene_name in self.TNGS_REGIONS.keys():
                self.logger.debug("LAB:iterate_section:[tNGS only] checking to see if this is a split primer")
                if isinstance(self.TNGS_REGIONS[variant.gene_name], dict):
                    for segment in self.TNGS_REGIONS[variant.gene_name]:
                        self.logger.debug("LAB:iterate_section:[tNGS only] checking if variant from {} is found in segment {}".format(variant.gene_name, segment))
                        if (self.TNGS_REGIONS[variant.gene_name][segment][0] <= variant.pos <= self.TNGS_REGIONS[variant.gene_name][segment][1]):
                            variant.gene_name_segment = segment
                            self.logger.debug("LAB:iterate_section:[tNGS only] variant from {} is found in segment {}; setting gene_name_segment to segment name".format(variant.gene_name, variant.gene_name_segment))
                            break

                    if hasattr(variant, "gene_name_segment") is False:
                        self.logger.warning("LAB:iterate_section:[tNGS only] This mutation is not in an expected region and is not within any segments")
                        variant.gene_name_segment = "Outside of expected region"

            # extract all of the annotations for the variant
            variant.extract_annotations()
            
            if variant.gene_name in ["mmpS5", "mmpL5", "mmpR5"]:
                reported_genes = variant.extract_consequences()            
                self.genes_reported.update(reported_genes)
            
            for annotation_row in variant.annotation_dictionary.values():
                # determine interpretations without QC overwrites    
                annotation_row.determine_interpretation()
            
                # perform QC on the row
                genes_with_valid_deletions, positional_qc_fails = annotation_row.add_qc_warnings(self.MIN_DEPTH, self.MIN_FREQUENCY, self.MIN_READ_SUPPORT, self.LOW_DEPTH_OF_COVERAGE_LIST, self.genes_with_valid_deletions)
                self.genes_with_valid_deletions.update(genes_with_valid_deletions)
                self.positional_qc_fails.update(positional_qc_fails)                
            
                # if in --debug mode, print the annotation row -- caution, this adds a ton of extra lines to the log.
                annotation_row.print()

                row_list.append(annotation_row)

        return row_list

    def create_laboratorian_report(self) -> pd.DataFrame:
        """
        This function creates the laboratorian report, which is a CSV file
        containing the following information for each mutation in the input JSON:
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
            - read_support: the number of reads supporting the mutation (10, depth*frequency)
            - rationale: the rationale for resistance calling (WHO classification, Expert rule)
            - warning: a column reserved for warnings such as low depth of coverage 
            - source: a column used to indicate the resistance source as specified by TBDB
            - tbdb_comment: a column used to include any additional comments as specified by TBDB
        """
        DF_LABORATORIAN = pd.DataFrame(columns = [
            "sample_id", "tbprofiler_gene_name", "tbprofiler_locus_tag", 
            "tbprofiler_variant_substitution_type", "tbprofiler_variant_substitution_nt",
            "tbprofiler_variant_substitution_aa", "confidence", "antimicrobial",
            "looker_interpretation", "mdl_interpretation", "depth", "frequency", 
            "read_support", "rationale", "warning", "source", "tbdb_comment"
        ])
        
        row_list = []
        self.logger.debug("LAB:create_laboratorian_report:Initializing the row_list; contains {} rows".format(len(row_list)))

        with open(self.input_json) as json_fh:
            input_json = json.load(json_fh)

            self.SAMPLE_NAME = input_json["id"]

            row_list = self.iterate_section(input_json["dr_variants"], row_list)
            row_list = self.iterate_section(input_json["other_variants"], row_list)

            self.logger.info("LAB:create_laboratorian_report:Iteration complete, there are now {} rows".format(len(row_list)))

        self.logger.debug("LAB:create_laboratorian_report:Now adding any genes that are missing from the report and editing/reordering any rows that need to be edited")
        for gene, antimicrobial_drug_names in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
            for drug_name in antimicrobial_drug_names:
                if gene not in self.genes_reported:
                    temporary_row = Row.wildtype_row(self.logger, self.SAMPLE_NAME, gene, drug_name, self.GENE_TO_LOCUS_TAG, self.GENE_TO_TIER, self.LOW_DEPTH_OF_COVERAGE_LIST)
                    row_list.append(temporary_row)

            # make a list to add QC fail rows to end of laboratorian report
            reorder_list = [] 

            if gene in self.LOW_DEPTH_OF_COVERAGE_LIST:
                for row in row_list:       
                    if row.tbprofiler_gene_name == gene:
                        if row.tbprofiler_gene_name in self.genes_with_valid_deletions:
                            if "Insufficient coverage in locus" in row.warning:
                                row.warning.remove("Insufficient coverage in locus")

                        # if the mutation fails quality control, move it to the end of the report
                        if ("Insufficient coverage in locus" in row.warning or row.tbprofiler_variant_substitution_nt in self.positional_qc_fails):
                            reorder_list.append(row)

                # remove rows in reorder_list from row_list and add them to the end of row_list
                for row in reorder_list:
                    row_list.remove(row)
                    row_list.append(row)

        self.logger.info("LAB:create_laboratorian_report:There are now {} rows after adding missing gene-drug combos and editing/reordering necessary rows\n".format(len(row_list)))
        
        for row in row_list:
            row.warning = list(filter(None, row.warning))
            row.warning = ". ".join(row.warning)

            # make a temporary dataframe out of the Row object using vars(row) which converts the object into a dictionary
            row_dictionary = pd.DataFrame(vars(row), index=[0], columns=DF_LABORATORIAN.columns)

            if len(DF_LABORATORIAN) == 0:
                DF_LABORATORIAN = row_dictionary
            else:
                DF_LABORATORIAN = pd.concat([DF_LABORATORIAN, row_dictionary], ignore_index=True)

        DF_LABORATORIAN.to_csv("{}.laboratorian_report.csv".format(self.output_prefix), index=False)
        
        return DF_LABORATORIAN