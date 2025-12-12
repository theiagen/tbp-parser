import globals as globals_
import pandas as pd
import json
import datetime

class LIMS:
    """
    This class creates the CDPH LIMS report.

    It has four functions:
        - get_id: returns the lineage in English for LIMS
        - convert_annotation: converts the resistance annotation and the target drug
            into the LIMS language
        - get_mutation_position: returns the position where a mutation occurs
        - apply_lims_rules: implements several parsing rules for the LIMS report
        - create_lims_report: creates the LIMS report CSV file
    """
    RESISTANCE_RANKING = {
        "R": 4,
        "U": 3,
        "S": 2,
        "WT": 1,
        "Insufficient Coverage": 0,
        "NA": -1 # outside the expected region
    }
    """A dictionary ranking the resistance annotations from highest to lowest priority."""
    
    SPECIAL_POSITIONS = {
        "rpoB": [426, 452], # codon; the RRDR range
        "gyrA": [88, 94], # codon; the QRDR range
        "gyrB": [446, 507], # codon; the QRDR range
        "rrl": [[2003, 2367], [2449, 3056]], # nucleotide; range
        "rrs": [1401, 1402, 1484] # nucleotide; specific positions
    }
    """This is a dictionary of positions for genes requiring different consideration.
    Note: the rpoB, gyrA, and gyrB special positions are in codons, rrl & rrs are 
    nucleotide positions, and all are ranges except rrs (indicates specific positions)
    """
    
    RPOB_MUTATIONS = [
        "Leu430Pro",
        "Asp435Tyr",
        "His445Asn",
        "His445Cys",
        "His445Leu",
        "His445Ser",
        "Leu452Pro",
        "Ile491Phe"
    ]
    """A list of rpoB mutations that require unique LIMS output wording."""

    def __init__(self, logger, input_json, output_prefix, LOW_DEPTH_OF_COVERAGE_LIST, SAMPLE_NAME, DF_LABORATORIAN, POSITIONAL_QC_FAILS, GENES_WITH_VALID_DELETIONS):
        self.logger = logger
        self.input_json = input_json
        self.output_prefix = output_prefix
                
        self.LOW_DEPTH_OF_COVERAGE_LIST = LOW_DEPTH_OF_COVERAGE_LIST
        self.SAMPLE_NAME = SAMPLE_NAME
        self.DF_LABORATORIAN = DF_LABORATORIAN
        self.POSITIONAL_QC_FAILS = POSITIONAL_QC_FAILS
        self.GENES_WITH_VALID_DELETIONS = GENES_WITH_VALID_DELETIONS

    def get_id(self, test=False) -> str:
        """Returns the lineage in English for LIMS

        Args:
            test (bool, optional): set to True if this is for a pytest. Defaults to False.

        Returns:
            str: the lineage in English
        """
        self.logger.info("LIMS:Within LIMS class get_id function")

        with open(self.input_json) as json_fh:
            input_json = json.load(json_fh)

            # set default values for lineage and sublineage from tbprofiler
            lineage = set()
            detected_lineage = input_json["main_lineage"]
            globals_.LINEAGE = detected_lineage
            detected_sublineage = input_json["sub_lineage"]
            self.logger.debug("LIMS:The detected lineage is: '{}', and the detected sublineage is: '{}'".format(detected_lineage, detected_sublineage))

            self.logger.debug("LIMS:Calculating the percentage of LIMS genes above the coverage threshold; the min locus percentage is set to {}".format(globals_.MIN_LOCUS_PERCENTAGE))

            try:
                if globals_.TNGS:
                    number_of_lims_genes_above_coverage_threshold = sum(float(globals_.COVERAGE_DICTIONARY[gene]) >= 90 for gene in globals_.COVERAGE_DICTIONARY.keys())
                    percentage_lims_genes_above = number_of_lims_genes_above_coverage_threshold / len(globals_.COVERAGE_DICTIONARY.keys())
                    self.logger.debug("LIMS:[tNGS only] The percentage of genes in the coverage dictionary above the coverage threshold is {};".format(percentage_lims_genes_above))
                else:
                    number_of_lims_genes_above_coverage_threshold = sum(float(globals_.COVERAGE_DICTIONARY[gene]) >= globals_.COVERAGE_THRESHOLD for gene in globals_.GENES_FOR_LIMS)
                    percentage_lims_genes_above = number_of_lims_genes_above_coverage_threshold / len(globals_.GENES_FOR_LIMS)

                self.logger.debug("LIMS:The percentage of LIMS genes above the coverage threshold is {}".format(percentage_lims_genes_above))

            except:
                self.logger.error("LIMS:Something went wrong -- this line shouldn't be printed unless a test is running! Setting percentage_lims_genes_above to 0")
                percentage_lims_genes_above = 0
                # lineage.add("DNA of Mycobacterium tuberculosis complex NOT detected")

            if test or (percentage_lims_genes_above >= globals_.MIN_LOCUS_PERCENTAGE):
                if globals_.TNGS:
                    self.logger.debug("LIMS:[tNGS only] The sequencing method is tNGS; now checking for a His57Asp mutation in pncA")
                    pncA_mutations = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == "pncA")]
                    if "p.His57Asp" in pncA_mutations["tbprofiler_variant_substitution_aa"].tolist():
                        self.logger.debug("LIMS:[tNGS only] p.His57Asp detected in pncA, lineage is likely M. bovis")
                        lineage.add("DNA of Mycobacterium bovis detected")
                    else:
                        self.logger.debug("LIMS:[tNGS only] p.His57Asp not detected in pncA, lineage is likely M. tuberculosis")
                        lineage.add("DNA of Mycobacterium tuberculosis complex detected (not M. bovis)") 

                else:  
                    self.logger.debug("LIMS:The sequencing method is WGS AND the percentage of LIMS genes is GREATER than {}% ({}); now checking the TBProfiler lineage calls".format(globals_.MIN_LOCUS_PERCENTAGE * 100, percentage_lims_genes_above))

                    sublineages = detected_sublineage.split(";")

                    if "lineage" in detected_lineage:
                        lineage.add("DNA of Mycobacterium tuberculosis species detected")

                    for sublineage in sublineages:  
                        if "BCG" in detected_lineage or "BCG" in sublineage:
                            lineage.add("DNA of Mycobacterium bovis BCG detected")

                        elif ("La1" in detected_lineage or "La1" in sublineage) or ("bovis" in detected_lineage or "bovis" in sublineage):
                            lineage.add("DNA of Mycobacterium bovis (not BCG) detected")     

                    if detected_lineage == "" or detected_lineage == "NA" or len(lineage) == 0:
                        self.logger.debug("LIMS:No lineage has been detected by TBProfiler; assuming M.tb")
                        lineage.add("DNA of Mycobacterium tuberculosis complex detected")

            else:
                self.logger.debug("LIMS:The percentage of LIMS genes is LESS than {}% ({}); assuming NOT M.tb".format(globals_.MIN_LOCUS_PERCENTAGE * 100, percentage_lims_genes_above))
                lineage.add("DNA of Mycobacterium tuberculosis complex NOT detected")

            lineage = "; ".join(sorted(lineage))

        globals_.LINEAGE_ENGLISH = lineage
        self.logger.debug("LIMS:The LIMS ID is: '{}'".format(globals_.LINEAGE_ENGLISH))
        self.logger.debug("LIMS:The TBProfiler lineage is: '{}'".format(globals_.LINEAGE))
        self.logger.info("LIMS:Finished getting lineage and ID, now exiting function")
        return lineage

    def convert_annotation(self, annotation, drug) -> str:
        """This function takes the resistance annotation and the target drug and converts it into the appropriate LIMS text

        Args:
            annotation (str): the highest resistance annotation associated with the drug
            drug (str): the drug for which the annotation is being converted

        Returns:
            str: the converted LIMS text for the drug
        """
        message = "No mutations associated with resistance to {} detected".format(drug)
        if annotation == "R":
            message = "Mutation(s) associated with resistance to {} detected".format(drug)
        elif (annotation == "R-Interim") or (annotation == "U"):
            message = "The detected mutation(s) have uncertain significance. Resistance to {} cannot be ruled out".format(drug)
        elif annotation == "Insufficient Coverage":
            message = "Pending Retest"      

        return message

    def apply_lims_rules(self, gene_dictionary, DF_LIMS, max_mdl_resistance, antimicrobial_code, responsible_gene, genes_associated_with_drug, drug, qc_pass_gene_drug_associated_rows) -> pd.DataFrame:
        """This function implements several parsing rules for the LIMS report. 

        Args:
            gene_dictionary (dict): the genes matched to their associated LIMS codes
            DF_LIMS (pd.DataFrame): the LIMS dataframe
            max_mdl_resistance (str): the maximum MDL resistance associated with this antimicrobial code
            antimicrobial_code (str): the antimicrobial code or header for the LIMS report
            responsible_gene (list[str]): a list of potential genes that may be responsible for the max_mdl_resistance
            genes_associated_with_drug (list[str]): a list of genes associated with the antimicrobial drug
            drug (str): the antimicrobial under question
            qc_pass_gene_drug_associated_rows (pd.DataFrame): the gene-drug associated rows that passed positional QC

        Returns:
            pd.DataFrame: the LIMS dataframe with the rules applied
        """ 
        self.logger.info("LIMS:Within LIMS class apply_lims_rules function; applying rules to mutations associated with {}".format(genes_associated_with_drug))

        mutations_per_gene = {}


        # i legitimately have no idea what's going on for the next 300 lines of code lol
         
        # get all MDL interpretations for the responsible gene(s) for the max mdl interpretations in case of recalculation
        all_responsible_mdl_interpretations = {}
        for responsible_gene_name in responsible_gene:
            all_responsible_mdl_interpretations[responsible_gene_name] = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == responsible_gene_name) & (globals_.DF_LABORATORIAN["antimicrobial"] == drug)]["mdl_interpretation"].tolist()

        self.logger.debug("LIMS:The MDL interpretations for the responsible gene(s) are: {}".format(all_responsible_mdl_interpretations))  

        for gene, gene_code in gene_dictionary.items():
            DF_LIMS[gene_code] = ""
            non_s_mutations = 0

            # check to see if we are reporting on this gene for the LIMS report
            gene_subset = qc_pass_gene_drug_associated_rows[qc_pass_gene_drug_associated_rows["tbprofiler_gene_name"] == gene]

            gene_subset = gene_subset[~gene_subset["warning"].str.contains("This mutation is outside the expected region")]

            # # get all mutations, their types, interpretations, and warnings associated with a particular gene and it's associated antimicrobial            
            # nt_mutations_per_gene = gene_subset["tbprofiler_variant_substitution_nt"].tolist()
            # aa_mutations_per_gene = gene_subset["tbprofiler_variant_substitution_aa"].tolist()
            # mutation_types_per_gene = gene_subset["tbprofiler_variant_substitution_type"].tolist()
            # mdl_interpretations = gene_subset["mdl_interpretation"].tolist()
            # warnings = gene_subset["warning"].tolist()
            # read_supports = gene_subset["read_support"].tolist()
            
            # Implementation of rule 5.3.3 - if there are any matching amino acid positions, keep the row with the higher read support
            self.logger.debug("LIMS:Considering if any mutations have identical amino acid positions and keeping only the one with higher read support")
            
            # make a subset so we don't have to deal with aggregating an entire dataframe after exploding
            duplicate_aa_check = gene_subset[["tbprofiler_variant_substitution_aa", "read_support"]].copy()
    
            # make a new column of the positions associated with each amino acid mutation and separate numeric and non-numeric positions
            duplicate_aa_check["aa_position"] = duplicate_aa_check["tbprofiler_variant_substitution_aa"].apply(lambda x: globals_.get_position(x))
            duplicate_aa_check_numeric = duplicate_aa_check[duplicate_aa_check["aa_position"].map(len) > 0].copy()
            duplicate_aa_check_non_numeric = duplicate_aa_check[duplicate_aa_check["aa_position"].map(len) == 0].copy()
            
            # explode the numeric positions so that each position has its own row
            duplicate_aa_check_exploded = duplicate_aa_check_numeric.explode("aa_position")
            
            # for each amino acid position, keep the row with the highest read support
            indices_to_keep = duplicate_aa_check_exploded.groupby("aa_position")["read_support"].idxmax()
            duplicate_aa_check_filtered = duplicate_aa_check_exploded.loc[indices_to_keep]
            
            # regroup by amino acid substitution to get back to original format
            duplicate_aa_check_filtered = duplicate_aa_check_filtered.groupby("tbprofiler_variant_substitution_aa", as_index=False).agg({"read_support": "first", "aa_positions": list})

            # add back in the non-numeric positions that were excluded earlier
            duplicate_aa_check_filtered = pd.concat([duplicate_aa_check_filtered, duplicate_aa_check_non_numeric], ignore_index=True)

            # keep only the rows in gene_subset that are in duplicate_aa_check_filtered -- we want to drop any that were filtered out
            rows_to_keep = gene_subset["tbprofiler_variant_substitution_aa"].isin(duplicate_aa_check_filtered["tbprofiler_variant_substitution_aa"])
            gene_subset_filtered = gene_subset[rows_to_keep].reset_index(drop=True)
            
            self.logger.debug("LIMS:After applying rule 5.3.3, the following amino acid mutations remain for this gene ({}) and drug ({}): {}".format(gene, drug, gene_subset_filtered["tbprofiler_variant_substitution_aa"].tolist()))
            
            for mutation in gene_subset_filtered["tbprofiler_variant_substitution_nt"].tolist():


                # grab the index so we can also grab the other information associated with this mutation
                index = nt_mutations_per_gene.index(mutation)
                aa_mutation = aa_mutations_per_gene[index]
                mutation_type = mutation_types_per_gene[index]

                # if it is rpoB, we want to get the position of the mutation to check if it is in RRDR
                if gene == "rpoB":           
                    position_aa = globals_.get_position(aa_mutation)

                # convert WT and Insufficient Coverage to empty strings
                if mutation == "WT" or mutation == "Insufficient Coverage" or mutation == "NA":
                    mutation = ""

                # convert NA,  WT, and Insufficient Coverage to empty strings
                if aa_mutation == "NA" or aa_mutation == "WT" or aa_mutation == "Insufficient Coverage":
                    aa_mutation = ""

                # do not add the mutation if the particular mutation has low quality or is blank
                if ("Failed quality in the mutation position" in warnings[index]) or ("Insufficient Coverage" in mdl_interpretations[index]) or (mutation == ""):
                    self.logger.debug("LIMS:This mutation (\"{}\", origin gene: {}) is not being added to the LIMS report because it failed quality in the mutation position, was WT, or had insufficient locus coverage".format(mutation, gene))
                    if "del" in mutation and "Failed quality in the mutation position" in warnings[index] and "Insufficient coverage in locus" in warnings[index]:
                        DF_LIMS[gene_code] = "No sequence"
                    else:
                        DF_LIMS[gene_code] = "No mutations detected"

                    # overwrite the MDL interpretation to be "WT" if the mutation is non-S
                    # see also rule 4.2.2.1
                    if self.RESISTANCE_RANKING[mdl_interpretations[index]] >= 2: 
                        self.logger.debug("LIMS:This mutation (\"{}\", origin gene: {}) is having its MDL interpretation rewritten to act as if WT".format(mutation, gene))
                        mdl_interpretations[index] = "WT"
                        if gene in responsible_gene:
                            all_responsible_mdl_interpretations[gene][index] = "WT"

                        if "del" in mutation and "Failed quality in the mutation position" in warnings[index]:
                            try:
                                if float(globals_.COVERAGE_DICTIONARY[gene]) < globals_.COVERAGE_THRESHOLD:
                                    mdl_interpretations[index] = "Insufficient Coverage"
                                    if gene in responsible_gene:
                                        all_responsible_mdl_interpretations[gene][index] = "Insufficient Coverage"
                                else:
                                    self.logger.debug("LIMS:This gene ({}) has sufficient coverage a deletion being present".format(gene))
                            except:
                                self.logger.debug("LIMS:This gene ({}) is not in the coverage dictionary".format(gene))

                        self.logger.debug("LIMS:Since this MDL interpretation changed, we are now potentially recalculating max_mdl_resistance (currently {})".format(max_mdl_resistance[0]))

                        if (len(all_responsible_mdl_interpretations) == 0):
                            self.logger.debug("LIMS:There are no more potential MDL interpretations; exiting loop")
                            break

                        if (max([self.RESISTANCE_RANKING[interpretation] for gene_set in all_responsible_mdl_interpretations.values() for interpretation in gene_set]) != self.RESISTANCE_RANKING[max_mdl_resistance[0]]) and gene in responsible_gene:
                            max_mdl_resistance = [annotation for annotation, rank in self.RESISTANCE_RANKING.items() if rank == max([self.RESISTANCE_RANKING[interpretation] for gene_set in all_responsible_mdl_interpretations.values() for interpretation in gene_set])]
                            self.logger.debug("LIMS:The maximum needed to be reevaluated; the potential new max_mdl_resistance is now {}".format(max_mdl_resistance[0]))

                            if DF_LIMS[antimicrobial_code][0] != "Pending Retest":
                                self.logger.debug("LIMS:Now changing the antimicrobial code value for this gene since max_mdl_resistance changed")
                                DF_LIMS[antimicrobial_code] = self.convert_annotation(max_mdl_resistance[0], drug)   
                            else:
                                self.logger.debug("LIMS:The antimicrobial code value is currently 'Pending Retest'; this shouldn't be overwritten.")
                        else:
                            self.logger.debug("LIMS:The maximum did not need to be reevaluated; the max_mdl_resistance is still {}".format(max_mdl_resistance[0]))
                # the mutation is of decent quality and non-S, we want to report all non-synonymous mutations UNLESS rpoB RRDR (see Variant l.145 for explanation)
                elif ((mutation_type != "synonymous_variant" and mdl_interpretations[index] != "S") 
                        or (gene == "rpoB" and globals_.is_within_range(position_aa, self.SPECIAL_POSITIONS[gene]))):
                    if aa_mutation != "" and aa_mutation != "p.0?": # report only amino acid mutations unless not applicable/blank/p.0?, in which case report nucleotide mutation
                        substitution = "{}".format(aa_mutation)
                    else:
                        substitution = "{}".format(mutation)

                    # the following if only captures synonymous mutations if rpoB RRDR mutations
                    if mutation_type == "synonymous_variant":
                        substitution = "{} [synonymous]".format(substitution)

                    # add the formatted mutation to mutations_per_gene dictionary (formatted as a list)
                    if gene not in mutations_per_gene.keys():
                        mutations_per_gene[gene] = [substitution]
                    elif substitution not in mutations_per_gene[gene]:
                        mutations_per_gene[gene] = "{}; {}".format("".join(mutations_per_gene[gene]), substitution)
                else:
                    self.logger.debug("LIMS:This mutation (\"{}\", origin gene: {}) is not being added to the LIMS report because it is not an rpoB RRDR \"S\" mutation".format(mutation, gene))

            # Mutations for a particular gene have been added to the mutations_per_gene dictionary.
            # if that gene has mutations associated with it, we want to perform some additional filtration,
            # especially if it is RRDR rpoB
            if gene in mutations_per_gene.keys():
                self.logger.debug("LIMS:Adding the mutations ({}) associated with this gene ({}) to the DF_LIMS dataframe".format(mutations_per_gene[gene], gene))
                DF_LIMS[gene_code] = mutations_per_gene[gene]

                if gene == "rpoB":
                    if max_mdl_resistance[0] == "R":
                        self.logger.debug("LIMS:This gene is rpoB, now checking to see if any of the mutations belong to the special position list")
                        rpob_specific_mutations_counter = 0

                        for mutation in mutations_per_gene[gene]:
                            if any(rpob_mutation in mutation for rpob_mutation in self.RPOB_MUTATIONS):
                                rpob_specific_mutations_counter += 1

                        if rpob_specific_mutations_counter > 0:
                            self.logger.debug("LIMS:Some rpoB mutations are in the special rpoB mutation list, changing the output message")
                            DF_LIMS[antimicrobial_code] = "Predicted low-level resistance to rifampin. May test susceptible by phenotypic methods."
                        else:
                            self.logger.debug("LIMS:None of the mutations are in the special rpoB mutation list, changing the output message")
                            DF_LIMS[antimicrobial_code] = "Predicted resistance to rifampin"

                    # if the *maximum* MDL resistance is S in this case and it is RRDR synonymous
                    elif max_mdl_resistance[0] == "S" and "synonymous" in mutations_per_gene[gene][0]:
                        self.logger.debug("LIMS:Only synonymous mutations were identified in RRDR; changing output message")
                        DF_LIMS[antimicrobial_code] = "Predicted susceptibility to rifampin. The detected synonymous mutation(s) do not confer resistance"

            elif gene == "rpoB":
                self.logger.debug("LIMS:No mutations were identified in rpoB; changing output message")
                DF_LIMS[antimicrobial_code] = "Predicted susceptibility to rifampin"

            else:
                self.logger.debug("LIMS:No mutations were identified for this gene-drug combination")

            # make sure that there is a mutation associated with this gene-drug combo
            if len(nt_mutations_per_gene) > 0:
                self.logger.debug("mdl_interpretations: {}".format(mdl_interpretations))
                maximum_ranking = max([self.RESISTANCE_RANKING[interpretation] for interpretation in mdl_interpretations])
                # see the RESISTANCE_RANKING dictionary for the ranking of each MDL interpretation
                # count up the number of non-S mutations; non-S mutations have a ranking > 2
                if maximum_ranking > 2:
                    non_s_mutations += 1

                # change the gene_code to be something different depending on the MDL interpretations and/or number of non-s mutations
                if maximum_ranking == 2 and non_s_mutations == 0 and ~ DF_LIMS[gene_code].str.contains("synonymous")[0]: # non-RRDR region S
                    DF_LIMS[gene_code] = "No high confidence mutations detected"
                elif maximum_ranking == 1: # WT
                    DF_LIMS[gene_code] = "No mutations detected"
                elif maximum_ranking == 0 and "del" not in DF_LIMS[gene_code]: # Insufficient Coverage
                    DF_LIMS[gene_code] = "No sequence"

            else:
                self.logger.debug("LIMS:There are no mutations for this gene ({}) associated with this drug ({})".format(gene, drug))
                DF_LIMS[gene_code] = "No mutations detected"

            if (DF_LIMS[gene_code][0] == "No sequence" or "Insufficient Coverage" in mdl_interpretations) and max_mdl_resistance[0] != "R" and gene not in globals_.GENES_WITH_DELETIONS:
                self.logger.debug("LIMS:This gene ({}) has insufficient coverage without a valid deletion and no other mutations associated with this antimicrobial ({}) are 'R'; changing antimicrobial output to 'Pending Retest'".format(gene, drug))
                DF_LIMS[antimicrobial_code] = "Pending Retest"    

        self.logger.info("LIMS:Finished applying rules to mutations associated with {}, exiting function".format(globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]))
        return DF_LIMS

    def create_lims_report(self) -> None:
        """
        This function recieves the input json file and laboratorian report to
        write the LIMS report that includes the following information:
            - MDL sample accession numbers:  sample name
            - M_DST_A01_ID - lineage
            - The set of information in ANTIMICROBIAL_CODE_TO_GENES dictionary with 
                target drug resistance information in layman's terms, and the 
                mutations responsible for the predicted phenotype
            - Date of analysis in YYYY-MM-DD HH:SS format
            - Operator information
        """
        self.logger.info("LIMS:Within LIMS class create_lims_report function")
        
        # TO-DO: make these column names configurable
        DF_LIMS = pd.DataFrame({
            "MDL sample accession numbers": self.SAMPLE_NAME, 
          "M_DST_A01_ID": self.get_id()
            }, index=[0])

        self.logger.debug("LIMS:Now iterating through each LIMS antimicrobial code")
        for drug, drug_gene_dictionary in globals_.ANTIMICROBIAL_CODE_TO_GENES.items():
            
            # get the MDL interpretations for all genes **FOR THE LIMS REPORT** associated with this drug       
            antimicrobial_code, gene_codes = drug_gene_dictionary.items()
            genes_associated_with_drug = gene_codes.keys()
            
            # capture all rows in DF_LABORATORIAN associated with this drug and then 
            #  capture only gene-drug associated rows that are reported on in the LIMS report
            drug_associated_rows = self.DF_LABORATORIAN[self.DF_LABORATORIAN["antimicrobial"] == drug]
            gene_drug_associated_rows = drug_associated_rows.loc[self.drug_associated_rows["tbprofiler_gene_name"].isin(genes_associated_with_drug)]          
                  
            # remove any interpretations that have failed positional QC from consideration
            # this is when the tbprofiler_variant_substitution_nt variable is in the POSITIONAL_QC_FAILS list
            qc_pass_gene_drug_associated_rows = gene_drug_associated_rows.loc[~gene_drug_associated_rows["tbprofiler_variant_substitution_nt"].isin(self.POSITIONAL_QC_FAILS)]
            
            potential_mdl_resistances = qc_pass_gene_drug_associated_rows["mdl_interpretation"].tolist()
            max_mdl_resistance = max(potential_mdl_resistances, key=lambda x: self.RESISTANCE_RANKING[x])
            
            self.logger.debug("LIMS:The max MDL resistance for this antimicrobial ({}) is {}".format(drug, max_mdl_resistance))

            # initalize list of genes responsible for the max resistance
            responsible_genes = set()

            try:
                # this will generate a list of annotations that match the value in maximum_resistance_identified
                rows_that_match_max = [self.gene_drug_associated_rows["mdl_interpretation"] == max_mdl_resistance]
                # get all gene names responsible for this max resistance
                responsible_genes = set(gene_drug_associated_rows.loc[rows_that_match_max]["tbprofiler_gene_name"].tolist())
                self.logger.debug("LIMS:The gene(s) responsible for the max MDL resistance for this antimicrobial ({}) is/are {}".format(drug, responsible_genes))
                
            except:
                max_mdl_resistance = "NA"
                
            DF_LIMS[antimicrobial_code] = self.convert_annotation(max_mdl_resistance, drug) 



            #### trying to remake apply_lim_rules here ####

            # don't apply any LIMS rules if the max MDL resistance is WT, Insufficient Coverage, or NA
            if max_mdl_resistance in ["WT", "Insufficient Coverage", "NA"]:
                # set all gene codes to "No mutations detected" or "No sequence" if Insufficient Coverage
                for gene, gene_code in gene_codes.items():
                    # if the gene has insufficient coverage without a valid deletion (NOT in GENES_WITH_VALID_DELETIONS but in LOW_DEPTH_OF_COVERAGE_LIST), set to "No sequence"
                    if gene in self.LOW_DEPTH_OF_COVERAGE_LIST and gene not in self.GENES_WITH_VALID_DELETIONS:
                        DF_LIMS[gene_code] = "No sequence"
                    else:
                        DF_LIMS[gene_code] = "No mutations detected"
                    
            if max_mdl_resistance == ["S"]:
                for gene, gene_code in gene_codes.items():
                    mutation_list = []
                    if gene == "rpoB" and drug == "rifampicin":
                        # check if there is a mutation appears in the RRDR (codons 426-452) (in the special-positions dictionary)
                        gene_subset = qc_pass_gene_drug_associated_rows[qc_pass_gene_drug_associated_rows["tbprofiler_gene_name"] == gene]
                        # reduce gene_subset to only those that are synonymous type
                        gene_subset = gene_subset[gene_subset["tbprofiler_variant_substitution_type"] == "synonymous_variant"]
                        if len(gene_subset) > 0:                    
                            for aa_mutation in gene_subset["tbprofiler_variant_substitution_aa"].tolist():
                                position_aa = globals_.get_position(aa_mutation)
                                
                                if globals_.is_within_range(position_aa, self.SPECIAL_POSITIONS[gene]):
                                    # synonymous RRDR mutation
                                    DF_LIMS[antimicrobial_code] = "Predicted susceptibility to rifampicin. The detected synonymous mutation(s) do not confer resistance"
                                    mutation_list.append("{} [synonymous]".format(aa_mutation))
                                else:
                                    # synonymous non-RRDR mutation -- ignore
                                    pass
                        else:
                            # no synonymous RRDR mutation but max MDL resistance is S and drug is rimpaciin
                            DF_LIMS[antimicrobial_code] = "Predicted susceptibility to rifampicin"
                        
                        DF_LIMS[gene_code] = "; ".join(mutation_list) if len(mutation_list) > 0 else "No high confidence mutations detected"
                    
                    else:
                        DF_LIMS[gene_code] = "No high confidence mutations detected"
            
            elif max_mdl_resistance in ["R", "U"]:
                # only report R or U mutations or RRDR synonymous mutations that are S
                for gene, gene_code in gene_codes.items():
                    mutation_list = []
                    gene_subset = qc_pass_gene_drug_associated_rows[qc_pass_gene_drug_associated_rows["tbprofiler_gene_name"] == gene]
                    
                    for index, row in gene_subset.iterrows():
                        mutation = row["tbprofiler_variant_substitution_nt"]
                        aa_mutation = row["tbprofiler_variant_substitution_aa"]
                        mutation_type = row["tbprofiler_variant_substitution_type"]
                        mdl_interpretation = row["mdl_interpretation"]
                        
                        if ((mdl_interpretation in ["R", "U"]) 
                            or (gene == "rpoB" and globals_.is_within_range(globals_.get_position(aa_mutation), self.SPECIAL_POSITIONS[gene]) and mutation_type == "synonymous_variant" and mdl_interpretation == "S")):
                            if aa_mutation != "" and aa_mutation != "p.0?": # report only amino acid mutations unless not applicable/blank/p.0?, in which case report nucleotide mutation
                                substitution = "{}".format(aa_mutation)
                            else:
                                substitution = "{}".format(mutation)

                            # the following if only captures synonymous mutations if rpoB RRDR mutations
                            if mutation_type == "synonymous_variant":
                                substitution = "{} [synonymous]".format(substitution)

                        mutation_list.append(substitution)
                        
                    if gene == "rpoB": # check special rpob positions to determine if we need to change teh antimicrobial code message
                        if max_mdl_resistance == "R":
                            rpob_specific_mutations_counter = 0

                            for mutation in mutation_list:
                                if any(rpob_mutation in mutation for rpob_mutation in self.RPOB_MUTATIONS):
                                    rpob_specific_mutations_counter += 1

                            if rpob_specific_mutations_counter > 0:
                                DF_LIMS[antimicrobial_code] = "Predicted low-level resistance to rifampicin. May test susceptible by phenotypic methods."
                            else:
                                DF_LIMS[antimicrobial_code] = "Predicted resistance to rifampicin"

                    if len(mutation_list) > 0:
                        DF_LIMS[gene_code] = "; ".join(mutation_list)
                    else:
                        DF_LIMS[gene_code] = "No mutations detected"
                            
                            
                        
                            
                                
                            
            
            # only apply LIMS rules if max MDL resistance is S for rpoB
            # apply LIMS rules if max MDL resistance is R or U
            
            

            # this function feels WAY too complicated and convoluted and i want to see if i can just remake it entirely
            #  especially the qc stuff and the rpob stuff -- i feel like we already capture a lot of the QC stuff elsewhere 
            #  and we're just repeating it needlessly here
            DF_LIMS = self.apply_lims_rules(gene_codes, DF_LIMS, max_mdl_resistance, antimicrobial_code, responsible_genes, genes_associated_with_drug, drug, qc_pass_gene_drug_associated_rows)

        DF_LIMS["Analysis date"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        DF_LIMS["Operator"] = globals_.OPERATOR

        # add lineage
        DF_LIMS["M_DST_O01_Lineage"] = globals_.LINEAGE

        # write to file
        DF_LIMS.to_csv("{}.lims_report.csv".format(self.output_prefix), index=False)
        self.logger.info("LIMS:LIMS report created, now exiting function\n")