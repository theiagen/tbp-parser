from Row import Row
from Variant import Variant
import globals as globals_
import pandas as pd
import json
import copy

class Laboratorian:
    """A class representing the Laboratorian report generation process.
    
    Attributes:
        logger (logging.Logger): The logger object for logging messages.
        input_json (str): The path to the input JSON file produced by TBProfiler.
        OUTPUT_PREFIX (str): The prefix for the output files.
        
        MIN_DEPTH (int): The minimum depth of coverage for QC.
        MIN_FREQUENCY (float): The minimum frequency for QC.
        MIN_READ_SUPPORT (int): The minimum read support for QC.
        COVERAGE_DICTIONARY (dict[str, float]): A dictionary containing coverage information for each gene/locus.
        LOW_DEPTH_OF_COVERAGE_LIST (list[str]): A list of genes/loci that have low depth of coverage.
        GENE_TO_ANTIMICROBIAL_DRUG_NAME (dict[str, str]): A dictionary mapping genes to their associated antimicrobial drug names.
        GENE_TO_LOCUS_TAG (dict[str, str]): A dictionary mapping genes to their locus tags.
        TNGS_REGIONS (dict[str, list[int] | dict[str, list[int]]]): A dictionary containing the specific tNGS regions for each gene (to account for non-overlapping primers in the same gene).
        GENE_TO_TIER (dict[str, str]): A dictionary mapping genes to their tiers.
        PROMOTER_REGIONS (dict[str, list[int] | list[list[int]]]): A dictionary containing promoter regions for each gene.
        TNGS (bool): A boolean indicating whether tNGS mode is enabled.
        
        genes_reported (set[str]): A set of genes that have already been reported in the laboratorian report.
        genes_with_valid_deletions (set[str]): A set of genes that have valid deletions reported in the laboratorian report.
        positional_qc_fails (dict[str, set[str]]): A dictionary of genes to mutations that have failed POSITIONAL QC checks.
    
    Methods:
        iterate_section(variant_section: dict, row_list: list[Row], raw_row_list: list[Row]) -> list[Row]:
            Iterates through each subsection in the "dr_variants" and "other_variants" sections
            in the TBProfiler JSON file, converting each item into a Variant object and extracting
            annotations into Row objects for the Laboratorian report.
        
        create_laboratorian_report() -> pd.DataFrame:
            Creates the laboratorian report as a CSV file containing mutation information
            extracted from the input JSON file.
    """

    def __init__(self, 
                 logger, 
                 input_json: str, 
                 OUTPUT_PREFIX: str, 
                 MIN_DEPTH: int, 
                 MIN_FREQUENCY: float, MIN_READ_SUPPORT: int, 
                 COVERAGE_DICTIONARY: dict[str, float], 
                 LOW_DEPTH_OF_COVERAGE_LIST: list[str], 
                 GENE_TO_ANTIMICROBIAL_DRUG_NAME: dict[str, str], 
                 GENE_TO_LOCUS_TAG: dict[str, str], 
                 TNGS_REGIONS: dict[str, list[int] | dict[str, list[int]]], 
                 GENE_TO_TIER: dict[str, str],
                 PROMOTER_REGIONS: dict[str, list[int] | list[list[int]]], 
                 TNGS: bool, DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY: bool) -> None:
        """Initializes the Laboratorian class with the provided parameters.
        
        Args:
            logger (logging.Logger): The logger object for logging messages.
            input_json (str): The path to the input JSON file produced by TBProfiler.
            OUTPUT_PREFIX (str): The prefix for the output files.
            MIN_DEPTH (int): The minimum depth of coverage for QC.
            MIN_FREQUENCY (float): The minimum frequency for QC.
            MIN_READ_SUPPORT (int): The minimum read support for QC.
            COVERAGE_DICTIONARY (dict[str, float]): A dictionary containing coverage information for each gene/locus.
            LOW_DEPTH_OF_COVERAGE_LIST (list[str]): A list of genes/loci that have low depth of coverage.
            GENE_TO_ANTIMICROBIAL_DRUG_NAME (dict[str, str]): A dictionary mapping genes to their associated antimicrobial drug names.
            GENE_TO_LOCUS_TAG (dict[str, str]): A dictionary mapping genes to their locus tags.
            TNGS_REGIONS (dict[str, list[int] | dict[str, list[int]]]): A dictionary containing the specific tNGS regions for each gene (to account for non-overlapping primers in the same gene).
            GENE_TO_TIER (dict[str, str]): A dictionary mapping genes to their tiers.
            PROMOTER_REGIONS (dict[str, list[int] | list[list[int]]]): A dictionary containing promoter regions for each gene.
            TNGS (bool): A boolean indicating whether tNGS mode is enabled.
            DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY (bool): A boolean indicating whether R mutations should be treated the same as S/U mutations for locus QC.
        """
        self.logger = logger
        self.input_json = input_json
        self.OUTPUT_PREFIX = OUTPUT_PREFIX
        
        # quality control variables
        self.MIN_DEPTH = MIN_DEPTH
        """An integer representing the minimum depth of coverage required for QC."""
        self.MIN_FREQUENCY = MIN_FREQUENCY  
        """A float representing the minimum frequency required for QC."""
        self.MIN_READ_SUPPORT = MIN_READ_SUPPORT
        """An integer representing the minimum read support required for QC."""
        
        # look-up dictionaries and lists 
        self.COVERAGE_DICTIONARY = COVERAGE_DICTIONARY
        """A dictionary containing coverage information for each gene/locus."""
        self.LOW_DEPTH_OF_COVERAGE_LIST = LOW_DEPTH_OF_COVERAGE_LIST
        """A list of genes/loci that have low depth of coverage."""
        self.GENE_TO_ANTIMICROBIAL_DRUG_NAME = GENE_TO_ANTIMICROBIAL_DRUG_NAME
        """A dictionary mapping genes to their associated antimicrobial drug names."""
        self.GENE_TO_LOCUS_TAG = GENE_TO_LOCUS_TAG
        """A dictionary mapping genes to their locus tags."""
        self.TNGS_REGIONS = TNGS_REGIONS
        """A dictionary containing the specific tNGS regions for each gene (to account for non-overlapping primers in the same gene)."""
        self.GENE_TO_TIER = GENE_TO_TIER
        """A dictionary mapping genes to their tiers."""
        self.PROMOTER_REGIONS = PROMOTER_REGIONS
        """A dictionary containing promoter regions for each gene."""
        self.TNGS = TNGS
        """A boolean indicating whether tNGS mode is enabled."""
        self.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY = DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY
        """A boolean indicating whether R mutations should be treated the same as S/U mutations for locus QC."""

        self.genes_reported = set()
        """A set of genes that have already been reported in the laboratorian report."""
        self.genes_with_valid_deletions = set()
        """A set of genes that have valid deletions reported in the laboratorian report."""
        self.positional_qc_fails = {}
        """A dictionary of genes to their mutations that have failed POSITIONAL QC checks."""

    def iterate_section(self, variant_section: dict, row_list: list[Row], raw_row_list: list[Row]) -> list[Row]:
        """This function iterates through each subsection in the "dr_variants" and "other_variants" sections
        in the TBProfiler JSON file. Each item in the subsection is converted into an individual Variant object, 
        and each annotation within that Variant is extracted and converted into a Row object, where each Row is
        a row in the Laboratorian report.

        Args:
            variant_section (json.load() object): The subsection of the TBProfiler JSON file to iterate through
            row_list (List): A list of rows that have been extracted so far
            raw_row_list (List): A list of rows (without QC applied) that have been extracted so far

        Returns:
            list[Row]: The list of rows, including the newly extracted rows
        """        
        self.logger.debug("LAB:iterate_section:Iterating through the variant section to turn each one into a Variant object")
        for variant in variant_section:
            # create a Variant object and add the origin gene to the class genes_reported set variable
            variant = Variant(self.logger, self.SAMPLE_NAME, 
                              self.GENE_TO_LOCUS_TAG, self.GENE_TO_TIER, 
                              self.GENE_TO_ANTIMICROBIAL_DRUG_NAME, 
                              self.PROMOTER_REGIONS, self.TNGS, variant)
            self.genes_reported.add(variant.gene_name)

            # extract all of the annotations for the variant
            variant.extract_annotations()
            
            variant.consequence_dictionary = []
            if variant.gene_name in ["mmpS5", "mmpL5", "mmpR5"]:
                reported_genes, annotation_dictionary = variant.extract_consequences()
                self.genes_reported.update(reported_genes)
                    
                variant.consequence_dictionary = annotation_dictionary

            combined_annotations = [variant.annotation_dictionary] + (variant.consequence_dictionary or [])
            for annotation in combined_annotations:
                for annotation_row in annotation.values():
                    # determine interpretations without QC overwrites    
                    annotation_row.determine_interpretation()
                
                    # make a copy of the annotation row that does not have any QC applied
                    raw_row = copy.deepcopy(annotation_row)
                    
                    # perform QC on the row
                    genes_with_valid_deletions, positional_qc_fails = annotation_row.add_qc_warnings(self.MIN_DEPTH, self.MIN_FREQUENCY, self.MIN_READ_SUPPORT, self.LOW_DEPTH_OF_COVERAGE_LIST, self.genes_with_valid_deletions, self.TNGS_REGIONS, self.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY)
                    self.genes_with_valid_deletions.update(genes_with_valid_deletions)
                    
                    if len(positional_qc_fails) > 0:
                        if variant.gene_name not in self.positional_qc_fails.keys():
                            self.positional_qc_fails[variant.gene_name] = positional_qc_fails
                        else:
                            self.positional_qc_fails[variant.gene_name].update(positional_qc_fails)

                    # if in --debug mode, print the annotation row -- caution, this adds a ton of extra lines to the log.
                    annotation_row.print()

                    raw_row_list.append(raw_row)
                    row_list.append(annotation_row)
        
        return row_list, raw_row_list
    
    def create_laboratorian_report(self) -> pd.DataFrame:
        """This function creates the laboratorian report, which is a CSV file
        containing the following information for each mutation in the input JSON:
            - sample_id: the sample name
            - tbprofiler_gene_name: the gene name
            - tbprofiler_locus_tag: the locus tag
            - tbprofiler_variant_substitution_type: the variant substitution type (missense_variant, upstream_gene_variant...)
            - tbprofiler_variant_substitution_nt: the nucleotide substitution (c.1349C>G)
            - tbprofiler_variant_substitution_aa: the amino acid substitution (p.Ser450Trp)
            - tbprofiler_variant_position: the position of the mutation in reference to the whole genome (761155)
            - confidence: the tbprofiler annotation regarding resistance (Not assoc w R, Uncertain significance...)
            - antimicrobial: the antimicrobial drug the mutation confers resistance to (streptomycin, rifampin...)
            - looker_interpretation: the interpretation of resistance for the CDPH Looker report (R, R-interim, U, S, S-interim)
            - mdl_interpretation: the MDL interpretation of resistance (R, S, U)
            - depth: the depth of coverage at the mutation site (100)
            - frequency: the frequency of mutation at the site (1)
            - read_support: the number of reads supporting the mutation (10, depth*frequency)
            - rationale: the rationale for resistance calling (WHO classification, Expert rule)
            - warning: a column reserved for warnings such as low depth of coverage 
            - gene_tier: the tier of the gene
            - source: a column used to indicate the resistance source as specified by TBDB
            - tbdb_comment: a column used to include any additional comments as specified by TBDB
            
        Returns:
            pd.DataFrame: The laboratorian report as a pandas DataFrame.
        """
        DF_LABORATORIAN = pd.DataFrame(columns = [
            "sample_id", "tbprofiler_gene_name", "tbprofiler_locus_tag", 
            "tbprofiler_variant_substitution_type", "tbprofiler_variant_substitution_nt",
            "tbprofiler_variant_substitution_aa", "tbprofiler_variant_position", "confidence", "antimicrobial",
            "looker_interpretation", "mdl_interpretation", "depth", "frequency", 
            "read_support", "rationale", "warning", "gene_tier", "source", "tbdb_comment"
        ])
        
        row_list = []
        raw_row_list = []
        self.logger.debug("LAB:create_laboratorian_report:Initializing the row_list; contains {} rows".format(len(row_list)))

        with open(self.input_json) as json_fh:
            input_json = json.load(json_fh)

            self.SAMPLE_NAME = input_json["id"]

            row_list, raw_row_list = self.iterate_section(input_json["dr_variants"], row_list, raw_row_list)
            row_list, raw_row_list = self.iterate_section(input_json["other_variants"], row_list, raw_row_list)

            self.logger.info("LAB:create_laboratorian_report:Iteration complete, there are now {} rows".format(len(row_list)))

        reorder_list = []
        self.logger.debug("LAB:create_laboratorian_report:Now adding any genes that are missing from the report and editing/reordering any rows that need to be edited")
        for gene, antimicrobial_drug_names in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
            
            for drug_name in antimicrobial_drug_names:
                if gene not in self.genes_reported:
                    temporary_row = Row.wildtype_row(self.logger, self.SAMPLE_NAME, gene, drug_name, self.GENE_TO_LOCUS_TAG, self.GENE_TO_TIER, self.LOW_DEPTH_OF_COVERAGE_LIST, self.TNGS)
                    row_list.append(temporary_row)
                    raw_row_list.append(copy.deepcopy(temporary_row))
                    
                    # add WT rows to bottom of report 
                    reorder_list.append(temporary_row)

            # add QC fail rows to end of laboratorian report too
            if gene in self.LOW_DEPTH_OF_COVERAGE_LIST:
                for row in row_list:       
                    if row.tbprofiler_gene_name == gene:
                        if row.tbprofiler_gene_name in self.genes_with_valid_deletions:
                            # remove insufficient coverage warning if a valid deletion was found
                            if "Insufficient coverage in locus" in row.warning:
                                row.warning.remove("Insufficient coverage in locus")

                        # if the mutation fails quality control, move it to the end of the report
                        if ("Insufficient coverage in locus" in row.warning 
                            or row.tbprofiler_variant_substitution_nt in self.positional_qc_fails.get(row.tbprofiler_gene_name, set())):
                            reorder_list.append(row)
                
        # reorder the rows so that WT and QC-fail rows are at the end of the report and sorted alphabetically by gene name
        remaining_rows = [row for row in row_list if row not in reorder_list]
        sorted_reorder_list = sorted(reorder_list, key=lambda x: x.tbprofiler_gene_name.lower())
        row_list = remaining_rows + sorted_reorder_list
        
        self.logger.info("LAB:create_laboratorian_report:There are now {} rows after adding missing gene-drug combos and editing/reordering necessary rows\n".format(len(row_list)))

        for rows in (row_list, raw_row_list):
            for row in rows:
                # concatenate warning list into a single string
                row.warning = list(filter(None, row.warning))
                row.warning = ". ".join(row.warning)

        DF_LABORATORIAN = pd.DataFrame([vars(row) for row in row_list], columns=DF_LABORATORIAN.columns)
        DF_LABORATORIAN.to_csv("{}.laboratorian_report.csv".format(self.OUTPUT_PREFIX), index=False)

        # the raw laboratorian report without any QC applied is here but will not be ordered the same as the main
        df_raw_laboratorian = pd.DataFrame([vars(row) for row in raw_row_list], columns=DF_LABORATORIAN.columns)
        df_raw_laboratorian.to_csv("{}.raw_lab_report.csv".format(self.OUTPUT_PREFIX), index=False)
                
        return DF_LABORATORIAN