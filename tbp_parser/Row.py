import globals as globals_

class Row() :
    """A class representing a row in the Laboratorian report.
    
    Attributes:
        sample_id (str): The sample ID.
        tbprofiler_gene_name (str): The gene name from TBProfiler.
        tbprofiler_locus_tag (str): The locus tag from TBProfiler.
        tbprofiler_variant_substitution_type (str): The type of variant substitution from TBProfiler.
        tbprofiler_variant_substitution_nt (str): The nucleotide change from TBProfiler.
        tbprofiler_variant_substitution_aa (str): The amino acid change from TBProfiler.
        confidence (str): The confidence used in the output report.
        who_confidence (str): The WHO confidence annotation.
        antimicrobial (str): The associated antimicrobial drug name.
        looker_interpretation (str): The Looker interpretation of the mutation.
        mdl_interpretation (str): The MDL interpretation of the mutation.
        depth (int): The depth of coverage at the mutation position.
        frequency (float): The allele frequency of the mutation.
        read_support (float): The read support for the mutation.
        rationale (str): The rationale for the interpretation.
        warning (set[str]): A set of warnings related to QC failures.
        gene_tier (str): The tier of the gene.
        source (str): The source of the annotation.
        tbdb_comment (str): Additional comments from the TBDB.
        variant (Variant): The Variant object associated with this row.
    
    Methods:
        wildtype_row(logger, sample_name: str, gene_name: str, drug_name: str, GENE_TO_LOCUS_TAG: dict[str, str], GENE_TO_TIER: dict[str, str], LOW_DEPTH_OF_COVERAGE_LIST: list[str]) -> 'Row':
            CLASS METHOD. Creates a wildtype Row object for genes that were not found in the TBProfiler
            JSON output (i.e., no mutations were found in that gene). This is used to create
            WT rows for genes that were sequenced but had no mutations.
            
        add_qc_warnings(MIN_DEPTH: int, MIN_FREQUENCY: float, MIN_READ_SUPPORT: float, LOW_DEPTH_OF_COVERAGE_LIST: list[str], genes_with_valid_deletions: set[str]) -> tuple[set, set]:
            Adds QC warnings if a mutation either has poor positional quality or locus quality.
            
        print() -> None:
            Prints the row in a readable format if debugging is enabled.
            
        determine_interpretation() -> None:
            Updates the Row's values for the rationale, confidence, and interpretations
            for the row.
            
        rank_annotation() -> int:
            Generates an integer ranking of the WHO annotation based on resistance,
            with 4 being the most resistant category and 1 the least.
        
        is_mutation_outside_region() -> bool:
            Checks if a mutation falls within the tNGS primer regions.
        
        describe_rationale(rule: str) -> str:
            Turns the rule number into a more descriptive string for the rationale
            column in the Laboratorian report.
    """
    
    def __init__(self, logger, variant: Variant, annotation: dict) -> None:
        """
        This function initializes the Row object with the appropriate
        values for each column in the CDPH Laboratorian report.
        No QC is performed.

        Args:
            logger: The logger object for logging messages.
            variant (Variant): A Variant object containing information about the mutation.
            annotation (dict): The annotation dictionary containing WHO confidence and drug information.
        """
        self.logger = logger

        self.variant = variant
        self.sample_id = variant.sample_name
        
        # extract fields from the variant object
        self.tbprofiler_gene_name = variant.__dict__.get("gene_name")

        self.depth = int(variant.__dict__.get("depth"))
        self.frequency = float(variant.__dict__.get("freq"))
        try:
            self.read_support = self.depth * self.frequency
        except:
            ### MATH FAILS
            raise Exception("MATH BAD")
        
        self.pos = variant.__dict__.get("pos")

        # should i use the .get() method here too?
        self.tbprofiler_variant_substitution_type = variant.type
        self.tbprofiler_variant_substitution_nt = variant.nucleotide_change
        self.tbprofiler_variant_substitution_aa = variant.protein_change if variant.protein_change != "" else "NA"

        # grab the locus_tag from the variant if it exists; otherwise, get it from the provided dictionary
        self.tbprofiler_locus_tag = variant.__dict__.get("locus_tag", variant.GENE_TO_LOCUS_TAG.get(self.tbprofiler_gene_name))

        self.gene_name = self.variant.__dict__.get("gene_name_segment", self.tbprofiler_gene_name)
        self.gene_tier = self.variant.GENE_TO_TIER.get(self.tbprofiler_gene_name, "NA")

        # extract fields from the annotation dictionary        
        self.who_confidence = annotation.get("confidence")
        if self.who_confidence == "" or annotation.get("comment") == "Not found in WHO catalogue":
            self.who_confidence = "No WHO annotation"
        self.confidence = self.who_confidence
        
        self.antimicrobial = annotation.get("drug")
        self.source = annotation.get("source", "")
        self.tbdb_comment = annotation.get("comment", "")

        # initialize empty values for the rest of the columns
        self.rationale = ""
        self.mdl_interpretation = ""
        self.looker_interpretation = ""
        self.warning = set()

    @classmethod
    def wildtype_row(cls, logger, sample_name: str, gene_name: str, drug_name: str, GENE_TO_LOCUS_TAG: dict[str, str], GENE_TO_TIER: dict[str, str], LOW_DEPTH_OF_COVERAGE_LIST: list[str]) -> 'Row':
        """This class method creates an wildtype Row object for genes that
        were not found in the TBProfiler JSON output (i.e., no mutations
        were found in that gene). This is used to create WT rows for
        genes that were sequenced but had no mutations.

        Args:
            logger: The logger object for logging messages.
            gene_name (str): The name of the gene.
            drug_name (str): The name of the associated drug
            GENE_TO_LOCUS_TAG (dict[str, str]): A dictionary mapping gene names to locus tags.
            GENE_TO_TIER (dict[str, str]): A dictionary mapping gene names to their tier.
            LOW_DEPTH_OF_COVERAGE_LIST (list[str]): A list of genes with low depth of coverage (failed locus QC).

        Returns:
            Row: An empty Row object with the appropriate NA or WT values.
        """        
        row = cls.__new__(cls)
        row.logger = logger # is this necessary? keeping it just in case
        row.variant = None
        row.sample_id = sample_name
        row.tbprofiler_gene_name = gene_name
        row.antimicrobial = drug_name
        row.who_confidence = "NA"
        row.confidence = "NA"
        row.depth = "NA"
        row.frequency = "NA"
        row.read_support = "NA"
        row.rationale = "NA"
        row.warning = set()
        row.source = ""
        row.tbdb_comment = ""
        row.tbprofiler_variant_substitution_type = "WT"
        row.tbprofiler_variant_substitution_nt = "WT"
        row.tbprofiler_variant_substitution_aa = "WT"
        row.tbprofiler_locus_tag = GENE_TO_LOCUS_TAG.get(gene_name, "NA")
        row.gene_tier = GENE_TO_TIER.get(gene_name, "NA")
        
        if gene_name in LOW_DEPTH_OF_COVERAGE_LIST: 
            # 4.2.2.3.1 - WT with insufficient coverage
            row.looker_interpretation = "Insufficient Coverage"
            row.mdl_interpretation = "Insufficient Coverage"
            row.warning.add("Insufficient coverage in locus")
        else:
            # 4.1 - WT with sufficient coverage
            row.looker_interpretation = "S"
            row.mdl_interpretation = "WT"
        
        return row

    def add_qc_warnings(self, MIN_DEPTH: int, MIN_FREQUENCY: float, MIN_READ_SUPPORT: float, LOW_DEPTH_OF_COVERAGE_LIST: list[str], genes_with_valid_deletions: set[str]) -> tuple[set[str], set[str]]:
        """Adds QC warnings if a mutation either has poor positional quality or locus quality

        Args:
            MIN_DEPTH (int): the minimum depth for a mutation to pass positional QC
            MIN_FREQUENCY (float): the minimum allele frequency for a mutation to pass positional QC
            MIN_READ_SUPPORT (float): the minimum read support for a mutation to pass positional QC 
            LOW_DEPTH_OF_COVERAGE_LIST (list[str]): a list of genes that have failed locus QC (i.e., breadth of coverage < min_percent_coverage)
            genes_with_valid_deletions (set[str]): a set of genes with deletions that pass positional QC

        Returns:
            tuple[set[str], set[str]]: the updated set of genes with valid deletions, and the set of mutations that failed positional QC
        """        
        positional_qc_fails = set()
        
        # checking positional qc now
        if (self.depth < MIN_DEPTH or self.frequency < MIN_FREQUENCY or self.read_support < MIN_READ_SUPPORT):
            if "del" not in self.tbprofiler_variant_substitution_nt: 
                # 4.2.1.1 - postiional qc fail; not a deletion
                positional_qc_fails.add(self.tbprofiler_variant_substitution_nt) 
                self.warning.add("Failed quality in the mutation position") 
                
            elif "del" in self.tbprofiler_variant_substitution_nt:
                if  (0 < self.depth and self.depth < MIN_DEPTH):
                    # 4.2.1.2 - postiional qc fail, deletion with some depth but not enough
                    positional_qc_fails.add(self.tbprofiler_variant_substitution_nt)
                    self.warning.add("Failed quality in the mutation position") 
                    
                elif (self.depth == 0 and self.frequency >= MIN_FREQUENCY): 
                    # 4.2.1.3 - deletion with zero depth but good frequency
                    pass
            
                elif (self.frequency < MIN_FREQUENCY):
                    # frequency is poor -- positional qc fail 
                    #### TO-DO: DO I NEED TO KEEP THIS? IT WAS DESCRIBED IN AN EMAIL BUT NOT IN THE INTERPRETATION DOCUMENT
                    positional_qc_fails.add(self.tbprofiler_variant_substitution_nt)
                    self.warning.add("Failed quality in the mutation position") 

        # checking locus qc now
        if self.gene_name in LOW_DEPTH_OF_COVERAGE_LIST:
            if "del" in self.tbprofiler_variant_substitution_nt: # 4.2.2.2 - locus qc fail and a deletion
                if self.tbprofiler_variant_substitution_nt in positional_qc_fails:
                    # this mutation also failed positional qc, so we need to add the locus warning
                    self.warning.add("Insufficient coverage in locus")
                    
            else: # 4.2.2.3 - locus qc fail but not a deletion -- if we maintain the "treat_r_as_s" option, we will implement that here
                if self.mdl_interpretation == "R" and "Failed quality in the mutation position" not in self.warning:
                    self.warning.add("Insufficient coverage in locus") # 4.2.2.3.3 - R mutation with locus qc fail but NOT positional qc fail; add warning DO NOT not overwrite interpretation
                
                elif self.mdl_interpretation == "R" and "Failed quality in the mutation position" in self.warning:
                    self.warning.add("Insufficient coverage in locus") # 4.2.2.3.4 - R mutation with BOTH positional and locus qc fail; add warning and overwrite interpretation
                    self.looker_interpretation = "Insufficient Coverage"
                    self.mdl_interpretation = "Insufficient Coverage"
                
                elif self.mdl_interpretation == "U" or self.mdl_interpretation == "S":
                    self.warning.add("Insufficient coverage in locus") # 4.2.2.3.2 - non-R mutation with locus qc fail; add warning and overwrite interpretation
                    self.looker_interpretation = "Insufficient Coverage"
                    self.mdl_interpretation = "Insufficient Coverage"
                
        if "Failed quality in the mutation position" not in self.warning and "del" in self.tbprofiler_variant_substitution_nt:
            genes_with_valid_deletions.add(self.gene_name)
            
        return genes_with_valid_deletions, positional_qc_fails

    def print(self) -> None:
        """This function prints the row in a readable format if debugging is enabled.
        """
        self.logger.debug("ROW:print:Now printing the row in a readable format:")
        self.logger.debug("ROW:print:\tsample_id: {}".format(self.sample_id))
        self.logger.debug("ROW:print:\ttbprofiler_gene_name: {}".format(self.tbprofiler_gene_name))
        self.logger.debug("ROW:print:\ttbprofiler_locus_tag: {}".format(self.tbprofiler_locus_tag))
        self.logger.debug("ROW:print:\ttbprofiler_variant_substitution_type: {}".format(self.tbprofiler_variant_substitution_type))
        self.logger.debug("ROW:print:\ttbprofiler_variant_substitution_nt: {}".format(self.tbprofiler_variant_substitution_nt))
        self.logger.debug("ROW:print:\ttbprofiler_variant_substitution_aa: {}".format(self.tbprofiler_variant_substitution_aa))
        self.logger.debug("ROW:print:\tconfidence: {}".format(self.confidence))
        self.logger.debug("ROW:print:\tantimicrobial: {}".format(self.antimicrobial))
        self.logger.debug("ROW:print:\tlooker_interpretation: {}".format(self.looker_interpretation))
        self.logger.debug("ROW:print:\tmdl_interpretation: {}".format(self.mdl_interpretation))
        self.logger.debug("ROW:print:\tdepth: {}".format(self.depth))
        self.logger.debug("ROW:print:\tfrequency: {}".format(self.frequency))
        self.logger.debug("ROW:print:\tread_support: {}".format(self.read_support))
        self.logger.debug("ROW:print:\trationale: {}".format(self.rationale))
        self.logger.debug("ROW:print:\twarning: {}".format(self.warning))
        self.logger.debug("ROW:print:\ttier: {}".format(self.gene_tier))
        self.logger.debug("ROW:print:\tsource: {}".format(self.source))
        self.logger.debug("ROW:print:\ttbdb_comment: {}".format(self.tbdb_comment))

    def determine_interpretation(self) -> None:
        """This function updates the Row's values for the rationale, confidence,
        and interpretations for the row. MDL interpretations are the same as Looker but
        drop the "- Interim" designations.
        """    
        ANNOTATION_TO_INTERPRETATION = {
            "Assoc w R": {
                "looker": "R", 
                "mdl": "R"
            },
            "Assoc w R - interim": {
                "looker": "R-Interim",
                "mdl": "R"
            },
            "Assoc w R - Interim": {
                "looker": "R-Interim",
                "mdl": "R"
            },
            "Uncertain significance": {
                "looker": "U", 
                "mdl": "U" 
            },
            "Not assoc w R": {
                "looker": "S",
                "mdl": "S"
            }, 
            "Not assoc w R - Interim": {
                "looker": "S-Interim", 
                "mdl": "S"
            }                              
        }
        
        if "This mutation is outside the expected region" in self.warning:
            self.rationale = "NA"
            self.confidence = "NA"
            self.looker_interpretation = "NA"
            self.mdl_interpretation = "NA"

        elif self.who_confidence != "No WHO annotation" and self.who_confidence != "" and self.who_confidence != "NA":
            self.looker_interpretation = ANNOTATION_TO_INTERPRETATION[self.who_confidence]["looker"]
            self.mdl_interpretation = ANNOTATION_TO_INTERPRETATION[self.who_confidence]["mdl"]
            self.rationale = "WHO classification"

        elif self.who_confidence != "NA":
            interpretation, rule = self.variant.apply_expert_rules()
            self.looker_interpretation = interpretation
            self.mdl_interpretation = interpretation
            self.rationale = self.describe_rationale(rule)
            self.confidence = "No WHO annotation"

    def rank_annotation(self) -> int: 
        """Generates an integer ranking of the WHO annotation based on resistance,
        with 4 being the most resistant category and 1 the least.
        
        Returns:
            int: The rank of the WHO annotation.
        """
        if self.who_confidence == "Assoc w R":
            return 4
        elif self.who_confidence == "Assoc w R - interim":
            return 3
        elif self.who_confidence == "Uncertain significance":
            return 2
        else:
            return 1  

    def is_mutation_outside_region(self) -> bool:
        """This function checks if a mutation falls within the primer regions.
            
            Returns:
                bool: True if the mutation is outside the expected region, false otherwise
        """
        # default to failure (mutation outside expected region)
        for primer, positions in globals_.TNGS_REGIONS.items():
            if primer != self.tbprofiler_gene_name:
                continue

            # split primers (positions is a dictionary)
            if isinstance(positions, dict):
                for segment, seg_positions in positions.items():
                    if segment != self.gene_name:
                        continue

                    if seg_positions[0] <= self.pos <= seg_positions[1]:
                        self.logger.debug("ROW:[tNGS only] Mutation falls within split primer region {}: {} is within {}".format(segment, self.pos, seg_positions))
                        return False  # position found, mutation is within region

                    self.logger.debug("ROW:[tNGS only] Mutation falls outside split primer region {}; {} is NOT within {}".format(segment, self.pos, seg_positions))

            # regular primers (positions is a list)
            else:
                if positions[0] <= self.pos <= positions[1]:
                    self.logger.debug("ROW:[tNGS only] Mutation falls within primer region {}: {} is within {}".format(primer, self.pos, positions))
                    return False  # position found, mutation is within region

                self.logger.debug("ROW:[tNGS only] Mutation does not fall within any primer regions for {}; {} is NOT within {}".format(primer, self.pos, positions))
                return True

        return True  # no matching position found, mutation is outside region

    def describe_rationale(self, rule: str) -> str:
        """This function turns the rule number into a more descriptive string for the rationale
        column in the Laboratorian report

        Args:
            rule (str): the rule as determined by Variant.apply_expert_rules()

        Returns:
            str: the rule with more description
        """        
        RULE_TO_RATIONALE = {
            "rule1.2": "Expert rule 1.2. Novel drug targets",
            "rule2.2.1": "Expert rule 2.2.1. Loss-of-function",
            "rule2.2.2.1": "Expert rule 2.2.2.1. rpoB RRDR",
            "rule2.2.2.2": "Expert rule 2.2.2.2 rpoB non-RRDR",
            "rule3.2.1": "Expert rule 3.2.1. rrs",
            "rule3.2.2": "Expert rule 3.2.3. gyrA QRDR",
            "rule3.2.3": "Expert rule 3.2.3. gyrB QRDR",
            "rule3.2.4": "No WHO annotation or expert rule",
            "whov2": "Mutation in proximal promoter region"
        }

        return RULE_TO_RATIONALE.get(rule, "No WHO annotation or expert rule")
      