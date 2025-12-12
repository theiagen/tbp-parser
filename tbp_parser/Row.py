import globals as globals_

class Row() :
    """
    This class represents a row in the Laboratorian report.
    
    The __init__ function assigns each variant attribute to the appropriate
    column in the Laboratorian report.

    This class has seven additional functions:
        - wildtype_row: creates a WT row for genes with no mutations
        - do_not_use: DEPRECATED - old init function that is no longer used
        - add_qc_warnings: adds QC warnings based on depth, frequency, and read 
           support
        - print: prints the row in a readable format
        - determine_interpretations: determines the Looker and MDL interpretations
        - rank_annotation: ranks the WHO annotation based on resistance
        - is_mutation_outside_region: checks if a mutation falls within the primer 
           regions
        - describe_rationale: describes the rationale for the interpretation if it
           originated from an expert rule
    """
    
    def __init__(self, logger, variant, annotation):
        """
        This function initializes the Row object with the appropriate
        values for each column in the CDPH Laboratorian report.
        No QC is performed.

        Args:
            logger: The logger object for logging messages.
            variant: A Variant object containing information about the mutation.
            annotation: The annotation dictionary containing WHO confidence and drug information.
        """
        self.logger = logger

        self.variant = variant
        self.sample_id = variant.sample_name
        
        # extract fields from the variant object
        self.tbprofiler_gene_name = variant.__dict__.get("gene_name")

        self.depth = int(variant.__dict__.get("depth"))
        self.frequency = float(variant.__dict__.get("freq"))
        try:
            self.read_support = int(self.depth * self.freq)
        except:
            ### MATH FAILS
            raise Exception("MATH BAD")
        
        self.pos = variant.__dict__.get("pos")

        # should i use the .get() method here too?
        self.tbprofiler_variant_substitution_type = variant.type
        self.tbprofiler_variant_substitution_nt = variant.nucleotide_change
        self.tbprofiler_variant_substitution_aa = variant.protein_change if variant.protein_change is not "" else "NA"

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
        self.warning = []

    @classmethod
    def wildtype_row(cls, logger, sample_name, gene_name, drug_name, GENE_TO_LOCUS_TAG, GENE_TO_TIER, LOW_DEPTH_OF_COVERAGE_LIST) -> 'Row':
        """This class method creates an wildtype Row object for genes that
        were not found in the TBProfiler JSON output (i.e., no mutations
        were found in that gene). This is used to create WT rows for
        genes that were sequenced but had no mutations.

        Args:
            logger: The logger object for logging messages.
            gene_name (str): The name of the gene.
            drug_name (str): The name of the associated drug.

        Returns:
            Row: An empty Row object with the appropriate NA or WT values.
        """        
        row = cls.__new1__(cls)
        row.logger = logger # is this necessary? keeping it just in case
        row.sample_name = sample_name
        row.tbprofiler_gene_name = gene_name
        row.antimicrobial = drug_name
        row.confidence = "NA"
        row.depth = "NA"
        row.frequency = "NA"
        row.read_support = "NA"
        row.rationale = "NA"
        row.warning = [""]
        row.source = ""
        row.tbdb_comment = ""
        row.tbprofiler_variant_substitution_type = "WT"
        row.tbprofiler_variant_substitution_nt = "WT"
        row.tbprofiler_variant_substitution_aa = "WT"
        row.locus_tag = GENE_TO_LOCUS_TAG.get(gene_name, "NA")
        row.gene_tier = GENE_TO_TIER.get(gene_name, "NA")
        
        if gene_name in LOW_DEPTH_OF_COVERAGE_LIST: 
            # 4.2.2.3.1 - WT with insufficient coverage
            row.looker_interpretation = "Insufficient Coverage"
            row.mdl_interpretation = "Insufficient Coverage"
            row.warning.append("Insufficient coverage in locus")
        else:
            # 4.1 - WT with sufficient coverage
            row.looker_interpretation = "S"
            row.mdl_interpretation = "WT"
        
        return row

    def do_not_use(self, logger, variant, who_confidence, drug, gene_name=None, depth=0, frequency=None, source="", tbdb_comment=""):
        self.logger = logger

        self.variant = variant
        self.who_confidence = who_confidence
        if tbdb_comment == "Not found in WHO catalogue":
            self.who_confidence = "No WHO annotation"
        self.antimicrobial = drug

        # create a row for the sample if the gene is in the coverage dictionary (should occur for all except in tNGS where only certain genes are sequenced)
        if gene_name != "test":
            self.antimicrobial = self.antimicrobial.replace("rifampicin", "rifampin")

            self.sample_id = globals_.SAMPLE_NAME

            # Initalizing the rest of the columns for the CDPH Laboratorian report
            # for when the variant is in the JSON file
            if variant is not None:
                self.logger.debug("ROW:__init__:Initalizing the Row object, the variant has been supplied.")

                try:
                    self.tbprofiler_gene_name = self.variant.gene_name
                except:
                    self.tbprofiler_gene_name = gene_name
                    self.variant.gene_name = gene_name
                if self.tbprofiler_gene_name == "mmpR5":
                    self.tbprofiler_gene_name = "Rv0678"

                self.logger.debug("ROW:__init__:The variant's gene name is {}".format(self.tbprofiler_gene_name))
                try:
                    self.tbprofiler_locus_tag = self.variant.locus_tag
                except:
                    self.tbprofiler_locus_tag = globals_.GENE_TO_LOCUS_TAG[self.tbprofiler_gene_name]
                self.tbprofiler_variant_substitution_type = self.variant.type
                self.tbprofiler_variant_substitution_nt = self.variant.nucleotide_change
                self.tbprofiler_variant_substitution_aa = self.variant.protein_change
                self.logger.debug("ROW:__init__:This mutation is a {} with nucleotide change \"{}\" and protein change \"{}\"".format(self.tbprofiler_variant_substitution_type, self.tbprofiler_variant_substitution_nt, self.tbprofiler_variant_substitution_aa))
                # change blank aa substitutions to NA
                if self.tbprofiler_variant_substitution_aa == "":
                    self.tbprofiler_variant_substitution_aa = "NA"
                self.confidence = self.who_confidence
                self.looker_interpretation = ""
                self.mdl_interpretation = ""
                try:
                    self.depth = int(self.variant.depth)
                except:
                    self.depth = depth
                try:
                    self.frequency = self.variant.freq
                except:
                    self.frequency = frequency
                # avoid division by zero errors
                try:
                    self.read_support = self.variant.depth * self.variant.freq
                except:
                    self.read_support = self.depth * self.frequency
                self.rationale = ""
                self.warning = []

                if hasattr(self.variant, "gene_name_segment"):
                    gene_name = self.variant.gene_name_segment
                else:
                    gene_name = self.tbprofiler_gene_name


                ### QC WILL BE SKIPPED FOR NOW ###
                if gene_name in globals_.COVERAGE_DICTIONARY:
                    if (self.depth < globals_.MIN_DEPTH) or (float(globals_.COVERAGE_DICTIONARY[gene_name]) < globals_.COVERAGE_THRESHOLD):
                        self.logger.debug("ROW:__init__:The depth of coverage for this variant is {} and the coverage for the gene is {}; applying a locus warning".format(self.depth, globals_.COVERAGE_DICTIONARY[gene_name]))
                        if (float(globals_.COVERAGE_DICTIONARY[gene_name]) < globals_.COVERAGE_THRESHOLD):
                            globals_.LOW_DEPTH_OF_COVERAGE_LIST.append(gene_name)

                            if "del" in self.tbprofiler_variant_substitution_nt or gene_name in globals_.GENES_WITH_DELETIONS:
                                self.logger.debug("ROW:__init__:This is a deletion, no warning added for the locus unless it fails positional qc (checked next)")
                            else:
                                self.warning.append("Insufficient coverage in locus")
                else:
                    self.logger.debug("ROW:__init__:This gene does not appear in the coverage dictionary. If tNGS, an additional warning will be given.")

                if globals_.TNGS:
                    self.logger.debug("ROW:__init__:[tNGS only] Checking to see if the mutation falls within the primer regions; if not, an additional warning will be given.")            
                    if self.is_mutation_outside_region(): # true if outside region
                        self.logger.debug("ROW:__init__:[tNGS only] This mutation's genomic position is outside the expected region")
                        self.warning.append("This mutation is outside the expected region")
                        self.logger.debug("ROW:__init__:[tNGS only] Rewriting this variant's interpretation to NA since it shouldn't exist")
                        self.looker_interpretation = "NA"
                        self.mdl_interpretation = "NA"

                protein_position = globals_.get_position(self.tbprofiler_variant_substitution_aa)

                # check to see if we need to apply a mutation warning 
                # (check rrs & rrl for low frequency and read support; 
                #  also check ethA & rpoB for specific protein position frequency)
                if ((self.depth < globals_.MIN_DEPTH) 
                     or ((self.tbprofiler_gene_name not in ["rrs", "rrl"] and 
                         (float(self.frequency) < globals_.MIN_FREQUENCY or self.read_support < globals_.MIN_READ_SUPPORT))
                     or (self.tbprofiler_gene_name == "rrs" and 
                         (float(self.frequency) < globals_.RRS_FREQUENCY or self.read_support < globals_.RRS_READ_SUPPORT)) 
                     or (self.tbprofiler_gene_name == "rrl" and 
                         (float(self.frequency) < globals_.RRL_FREQUENCY or self.read_support < globals_.RRL_READ_SUPPORT)) 
                     or (self.tbprofiler_gene_name == "ethA" and 
                         237 in protein_position and float(self.frequency) < globals_.ETHA237_FREQUENCY)
                     or (self.tbprofiler_gene_name == "rpoB" and 
                         449 in protein_position and float(self.frequency) < globals_.RPOB449_FREQUENCY))):

                    if ("del" in self.tbprofiler_variant_substitution_nt and self.depth == 0 and self.read_support == 0):
                        # placeholder
                        self.logger.debug("ROW:__init__:This is an okay scenario")
                    else:
                        self.logger.debug("ROW:__init__:The depth of coverage for this variant is {}, the frequency is {}, and the read support is {}; applying an additional mutation position warning".format(self.depth, self.frequency, self.read_support))

                        if self.tbprofiler_gene_name in globals_.COVERAGE_DICTIONARY.keys():
                            if ((float(globals_.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) < globals_.COVERAGE_THRESHOLD) and 
                                ("del" in self.tbprofiler_variant_substitution_nt or self.tbprofiler_gene_name in globals_.GENES_WITH_DELETIONS)):
                                self.logger.debug("ROW:__init__:This deletion failed in the mutation position and there was insufficient coverage locus, adding insufficient coverage warning")
                                self.warning.append("Insufficient coverage in locus")

                        globals_.MUTATION_FAIL_LIST.append(self.tbprofiler_variant_substitution_nt)
                        self.warning.append("Failed quality in the mutation position")

                elif globals_.TNGS:
                    self.logger.debug("ROW:__init__:[tNGS only] Checking boundaries for QC")
                    lower_rs = globals_.TNGS_READ_SUPPORT_BOUNDARIES[0]
                    upper_rs = globals_.TNGS_READ_SUPPORT_BOUNDARIES[1]
                    lower_f = globals_.TNGS_FREQUENCY_BOUNDARIES[0]
                    upper_f = globals_.TNGS_FREQUENCY_BOUNDARIES[1]

                    tngs_fail = False
                    if (lower_rs <= self.read_support and self.read_support < upper_rs):
                        if (self.frequency >= upper_f):
                            self.logger.debug("ROW:__init__:[tNGS only] this mutation's read support ({}) is between {} and {} and the frequency is above {}, so it passed the tNGS QC boundaries, no warning added for the mutation position".format(self.read_support, lower_rs, upper_rs, upper_f))
                        else:
                            self.logger.debug("ROW:__init__:[tNGS only] this mutation's read support ({}) is between {} and {} but the frequency is below {}, so it failed the tNGS QC boundaries".format(self.read_support, lower_rs, upper_rs, upper_f))
                            tngs_fail = True
                    elif (self.read_support >= upper_rs):
                        if (self.frequency >= lower_f):
                            self.logger.debug("ROW:__init__:[tNGS only] this mutation's read support ({}) is above {} and the frequency is above {}, so it passed the tNGS QC boundaries, no warning added for the mutation position".format(self.read_support, upper_rs, lower_f))
                        else:
                            self.logger.debug("ROW:__init__:[tNGS only] this mutation's read support ({}) is above {} but the frequency is below {}, so it failed the tNGS QC boundaries".format(self.read_support, upper_rs, lower_f))
                            tngs_fail = True
                    elif (self.read_support < lower_rs):
                        self.logger.debug("ROW:__init__:[tNGS only] this mutation's read support ({}) is below {}, so it failed the tNGS QC boundaries".format(self.read_support, lower_rs))
                        tngs_fail = True

                    if tngs_fail:
                        globals_.MUTATION_FAIL_LIST.append(self.tbprofiler_variant_substitution_nt)
                        self.warning.append("Failed quality in the mutation position")

                else:
                    self.logger.debug("ROW:__init__:The depth of coverage for this variant is {}, the frequency is {}, and the read support is {}; no additional warning added for the mutation position".format(self.depth, self.frequency, self.read_support))
                    if len(self.warning) == 0:
                        self.warning = [""]

                self.logger.debug("ROW:__init__:This variant has the following warnings: {}".format(self.warning))

                if "Failed quality in the mutation position" not in self.warning and "del" in self.tbprofiler_variant_substitution_nt:
                    globals_.GENES_WITH_DELETIONS.add(self.tbprofiler_gene_name)
                    self.logger.debug("ROW:__init__:This is a deletion that passed positional qc, adding to set")

            # otherwise, the variant does not appear in the JSON file (or was not sequenced [tNGS]) and default NA/WT values need to be supplied
            else:
                self.logger.debug("ROW:__init__:Initializing the Row object, the variant has no information supplied. Defaulting to NA or WT values.")

                if gene_name == "mmpR5":
                    self.tbprofiler_gene_name = "Rv0678"
                else: 
                    self.tbprofiler_gene_name = gene_name
                try:
                    self.tbprofiler_locus_tag = globals_.GENE_TO_LOCUS_TAG[self.tbprofiler_gene_name]
                except:
                    self.tbprofiler_locus_tag = "NA"
                self.confidence = "NA"
                self.depth = "NA"
                self.frequency = "NA"
                self.read_support = "NA"
                self.rationale = "NA"
                self.warning = [""]

                # check to see if we need to apply a coverage warning 
                #  (whole locus fail point c; rule 4.1 & 4.2.1.3.1)
                if (self.tbprofiler_gene_name in globals_.COVERAGE_DICTIONARY.keys() 
                    or self.tbprofiler_gene_name in globals_.TNGS_REGIONS.keys()):
                    try: 
                        if float(globals_.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) >= globals_.COVERAGE_THRESHOLD:
                            self.logger.debug("ROW:__init__:A warning does not need to be applied; setting the variant information to WT")
                            self.tbprofiler_variant_substitution_type = "WT"
                            self.tbprofiler_variant_substitution_nt = "WT"
                            self.tbprofiler_variant_substitution_aa = "WT"
                            self.looker_interpretation = "S"
                            self.mdl_interpretation = "WT"
                        else:
                            self.logger.debug("ROW:__init__:A warning needs to be applied; setting the variant information to insufficient coverage")
                            self.tbprofiler_variant_substitution_type = "Insufficient Coverage"
                            self.tbprofiler_variant_substitution_nt = "WT"
                            self.tbprofiler_variant_substitution_aa = "WT"
                            self.looker_interpretation = "Insufficient Coverage"
                            self.mdl_interpretation = "Insufficient Coverage"
                            self.warning.append("Insufficient coverage in locus")
                    except:
                        self.logger.debug("ROW:__init__:[tNGS only]: The gene does not appear in the coverage dictionary, but is in the TNGS regions dictionary")
                        self.logger.debug("ROW:__init__:[tNGS only]: This indicates that the gene was sequenced, but coverage was calculated under a different name")
                        self.logger.debug("ROW:__init__:[tNGS only]: A coverage warning will be applied if at least one segment is under the coverage threshold")

                        for segment in globals_.TNGS_REGIONS[self.tbprofiler_gene_name]:
                            if float(globals_.COVERAGE_DICTIONARY[segment]) >= globals_.COVERAGE_THRESHOLD:
                                self.logger.debug("ROW:__init__:[tNGS only]: This segment has good coverage, checking the other segments.")
                                self.tbprofiler_variant_substitution_type = "WT"
                                self.tbprofiler_variant_substitution_nt = "WT"
                                self.tbprofiler_variant_substitution_aa = "WT"
                                self.looker_interpretation = "S"
                                self.mdl_interpretation = "WT"
                            else:
                                self.logger.debug("ROW:__init__:[tNGS only]: This segment has poor coverage and a warning needs to be applied; setting the variant information to insufficient coverage")
                                self.tbprofiler_variant_substitution_type = "Insufficient Coverage"
                                self.tbprofiler_variant_substitution_nt = "WT"
                                self.tbprofiler_variant_substitution_aa = "WT"
                                self.looker_interpretation = "Insufficient Coverage"
                                self.mdl_interpretation = "Insufficient Coverage"
                                self.warning.append("Insufficient coverage in locus")
                                break

                else:
                    self.logger.debug("ROW:__init__:This gene ({}) was not sequenced".format(self.tbprofiler_gene_name))
                    self.tbprofiler_variant_substitution_type = "NA"
                    self.tbprofiler_variant_substitution_nt = "NA"
                    self.tbprofiler_variant_substitution_aa = "NA"
                    self.looker_interpretation = "NA"
                    self.mdl_interpretation = "NA"

            try:
                self.gene_tier = globals_.GENE_TO_TIER[self.tbprofiler_gene_name]
            except:
                try:
                    # iterate through the tNGS regions to see if we have a match
                    parent_genes = globals_.TNGS_REGIONS.keys()
                    parent_gene = [gene for gene in parent_genes if self.variant.gene_name_segment in globals_.TNGS_REGIONS[gene].keys()][0]
                    self.logger.debug("ROW:__init__:[tNGS only]: The parent gene ({}) of this segment ({}) was identified; now adding tier".format(parent_gene, self.variant.gene_name_segment))
                    self.gene_tier = globals_.GENE_TO_TIER[parent_gene]
                    # # now that coverage has been calculated, we can now rename the gene to be the parent gene name
                    # self.logger.debug("ROW:[tNGS only]: Renaming the gene segment ({}) to be the parent gene name ({})".format(self.tbprofiler_gene_name, parent_gene))
                    # self.tbprofiler_gene_name = parent_gene

                except:
                    self.gene_tier = "NA"

            self.source = source
            self.tbdb_comment = tbdb_comment

    def add_qc_warnings(self, MIN_DEPTH, MIN_FREQUENCY, MIN_READ_SUPPORT, LOW_DEPTH_OF_COVERAGE_LIST, genes_with_valid_deletions):
        positional_qc_fails = []
        
        # checking positional qc now
        if (self.depth < MIN_DEPTH or self.frequency < MIN_FREQUENCY or self.read_support < MIN_READ_SUPPORT):
            if "del" not in self.tbprofiler_variant_substitution_nt: 
                # 4.2.1.1 - postiional qc fail; not a deletion
                positional_qc_fails.append(self.tbprofiler_variant_substitution_nt) 
                self.warning.append("Failed quality in the mutation position") 
                
            elif "del" in self.tbprofiler_variant_substitution_nt:
                if  (0 < self.depth and self.depth < MIN_DEPTH):
                    # 4.2.1.2 - postiional qc fail, deletion with some depth but not enough
                    positional_qc_fails.append(self.tbprofiler_variant_substitution_nt)
                    self.warning.append("Failed quality in the mutation position") 
                    
                elif (self.depth == 0 and self.frequency >= MIN_FREQUENCY): 
                    # 4.2.1.3 - deletion with zero depth but good frequency
                    pass
            
                elif (self.frequency < MIN_FREQUENCY):
                    # frequency is poor -- positional qc fail 
                    #### TO-DO: DO I NEED TO KEEP THIS? IT WAS DESCRIBED IN AN EMAIL BUT NOT IN THE INTERPRETATION DOCUMENT
                    positional_qc_fails.append(self.tbprofiler_variant_substitution_nt)
                    self.warning.append("Failed quality in the mutation position") 

        # checking locus qc now
        if self.gene_name in LOW_DEPTH_OF_COVERAGE_LIST:
            if "del" in self.tbprofiler_variant_substitution_nt: # 4.2.2.2 - locus qc fail and a deletion
                if self.tbprofiler_variant_substitution_nt in positional_qc_fails:
                    # this mutation also failed positional qc, so we need to add the locus warning
                    self.warning.append("Insufficient coverage in locus")
                    
            else: # 4.2.2.3 - locus qc fail but not a deletion -- if we maintain the "treat_r_as_s" option, we will implement that here
                if self.mdl_interpretation == "R" and "Failed quality in the mutation position" not in self.warning:
                    self.warning.append("Insufficient coverage in locus") # 4.2.2.3.3 - R mutation with locus qc fail but NOT positional qc fail; add warning DO NOT not overwrite interpretation
                
                elif self.mdl_interpretation == "R" and "Failed quality in the mutation position" in self.warning:
                    self.warning.append("Insufficient coverage in locus") # 4.2.2.3.4 - R mutation with BOTH positional and locus qc fail; add warning and overwrite interpretation
                    self.looker_interpretation = "Insufficient Coverage"
                    self.mdl_interpretation = "Insufficient Coverage"
                
                elif self.mdl_interpretation == "U" or self.mdl_interpretation == "S":
                    self.warning.append("Insufficient coverage in locus") # 4.2.2.3.2 - non-R mutation with locus qc fail; add warning and overwrite interpretation
                    self.looker_interpretation = "Insufficient Coverage"
                    self.mdl_interpretation = "Insufficient Coverage"

                
        if "Failed quality in the mutation position" not in self.warning and "del" in self.tbprofiler_variant_substitution_nt:
            genes_with_valid_deletions.add(self.gene_name)
            
        return genes_with_valid_deletions, positional_qc_fails

    def print(self) -> None:
        """
        This function prints the row in a readable format.
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
        """
        This function finishes each row with the rest of the values needed. It depends
        on the ANNOTATION_TO_INTERPRETATION dictionary which is a dictionary to turn
        TBProfiler WHO annotations into their corresponding Looker or MDL
        interpretations; MDL interpretations are the same as Looker but drop the 
        "- Interim" designations.
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
            self.logger.debug("ROW:determine_interpretation:This mutation shouldn't exist! Setting Looker & MDL interpretations of 'NA'")
            self.rationale = "NA"
            self.confidence = "NA"
            self.looker_interpretation = "NA"
            self.mdl_interpretation = "NA"

        elif self.who_confidence != "No WHO annotation" and self.who_confidence != "" and self.who_confidence != "NA":
            self.logger.debug("ROW:determine_interpretation:WHO annotation identified: converting to the appropriate interpretation")
            self.looker_interpretation = ANNOTATION_TO_INTERPRETATION[self.who_confidence]["looker"]
            self.mdl_interpretation = ANNOTATION_TO_INTERPRETATION[self.who_confidence]["mdl"]
            self.rationale = "WHO classification"

        elif self.who_confidence != "NA":
            self.logger.debug("ROW:determine_interpretation:No WHO annotation identified: convert with expert rules")
            
            interpretation, rule = self.variant.apply_expert_rules()
            self.looker_interpretation = interpretation
            self.mdl_interpretation = interpretation
            self.rationale = self.describe_rationale(rule)
            self.confidence = "No WHO annotation"

    def rank_annotation(self) -> int: 
        """
        This function ranks the WHO annotation based on resistance,
        with 4 being the most resistant category and 1 the least.
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
        """
        This function checks if a mutation falls within the primer regions.
        Returns a boolean flag indicating if the mutation is outside the expected region.
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

    def describe_rationale(self, rule) -> str:
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
      