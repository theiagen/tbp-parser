import globals as globals_
import pandas as pd
import json
import datetime

class LIMS:
    """A class to perform LIMS-report specific functions
    
    Attributes:
        logger: the logger object
        input_json (str): the path to the TBProfiler JSON file
        OUTPUT_PREFIX (str): the prefix for output files
        LOW_DEPTH_OF_COVERAGE_LIST (list[str]): a list of genes with low depth of coverage (failed locus QC)
        SAMPLE_NAME (str): the sample name
        DF_LABORATORIAN (pd.DataFrame): a pandas DataFrame containing the laboratorian report data
        POSITIONAL_QC_FAILS (dict[str, str]): a dictionary mapping genes to their nucleotide mutation that failed positional QC
        GENES_WITH_VALID_DELETIONS (dict[str, list[int]]): a dictionary of genes (to the deletion genomic positions) with valid (QC-pass) deletions
        LINEAGE (str): the TBProfiler identified lineage
        LINEAGE_ENGLISH (str): the lineage in human-friendly language
        
        RESISTANCE_RANKING (dict[str, int]): a dictionary ranking the resistance annotations from highest to lowest priority
        SPECIAL_POSITIONS (dict[str, list[int] | list[list[int]]]): a dictionary of positions for genes requiring different consideration
        RPOB_MUTATIONS (list[str]): a list of rpoB mutations that require unique LIMS output wording
    
    Methods:
        get_id(TNGS: bool, MIN_PERCENT_LOCI_COVERED: float, test: bool = False) -> str:
            Returns both the TBProfiler identified lineage and the lineage in human-friendly language
            
        convert_annotation(annotation, drug) -> str:
            Converts the resistance annotation into the appropriate LIMS text
            
        create_lims_report(TNGS: bool, MIN_PERCENT_LOCI_COVERED: float, OPERATOR: str) -> None:
            Creates the LIMS report using the laboratorian pd.DataFrame
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

    def __init__(self, logger, input_json: str, OUTPUT_PREFIX: str, LOW_DEPTH_OF_COVERAGE_LIST: list[str], SAMPLE_NAME: str, DF_LABORATORIAN: pd.DataFrame, POSITIONAL_QC_FAILS: dict[str, str], GENES_WITH_VALID_DELETIONS: dict[str, list[int]]) -> None:
        """Initializes the LIMS class
        
        Args:
            logger: the logger object
            input_json (str): the path to the TBProfiler JSON file
            OUTPUT_PREFIX (str): the prefix for output files
            LOW_DEPTH_OF_COVERAGE_LIST (list[str]): a list of genes with low depth of coverage (failed locus QC)
            SAMPLE_NAME (str): the sample name
            DF_LABORATORIAN (pd.DataFrame): a pandas DataFrame containing the laboratorian report data
            POSITIONAL_QC_FAILS (dict[str, str]): a dictionary mapping genes to their nucleotide mutation that failed positional QC
            GENES_WITH_VALID_DELETIONS (dict[str, list[int]]): a dictionary of genes (to the deletion genomic positions) with valid (QC-pass) deletions
        """
        self.logger = logger
        self.input_json = input_json
        self.OUTPUT_PREFIX = OUTPUT_PREFIX

        self.LOW_DEPTH_OF_COVERAGE_LIST = LOW_DEPTH_OF_COVERAGE_LIST
        """A list of genes with low depth of coverage (failed locus QC)."""
        self.SAMPLE_NAME = SAMPLE_NAME
        """The sample name"""
        self.DF_LABORATORIAN = DF_LABORATORIAN
        """A pandas DataFrame containing the laboratorian report data."""
        self.POSITIONAL_QC_FAILS = POSITIONAL_QC_FAILS
        """A dictionary mapping genes to their nucleotide mutation that failed positional QC"""
        self.GENES_WITH_VALID_DELETIONS = GENES_WITH_VALID_DELETIONS
        """A dictionary of genes (to the deletion genomic positions) with valid (QC-pass) deletions"""

    def get_id(self, TNGS: bool, MIN_PERCENT_LOCI_COVERED: float, test: bool = False) -> str:
        """Returns both the TBProfiler identified lineage and the lineage in human-friendly languager

        Args:
            TNGS (bool): flag to indicate if the input data is tNGS or not
            MIN_PERCENT_LOCI_COVERED (float): the minimum percentage of the loci that need to be covered to call a lineage
            test (bool, optional): determines if this is running as part of a pytest module. Defaults to False.

        Returns:
            str: the lineage in English
        """
        with open(self.input_json) as json_fh:
            input_json = json.load(json_fh)

            # set default values for lineage and sublineage from tbprofiler
            lineage = set()
            detected_lineage = input_json["main_lineage"]
            self.LINEAGE = detected_lineage
            detected_sublineage = input_json["sub_lineage"]
                       
            lims_genes = list({
                gene_name for code_to_genes in globals_.DRUG_COLUMNS_TO_GENE_COLUMNS.values()
                for gene_to_code in code_to_genes.values()
                for gene_name in gene_to_code.keys()
            })
            
            # count number of LIMS genes in the LOW_DEPTH_OF_COVERAGE_LIST that do not have valid deletions that may explain the coverage
            passing_gene_count = 0
            for gene in lims_genes:
                if gene not in self.LOW_DEPTH_OF_COVERAGE_LIST and gene not in self.GENES_WITH_VALID_DELETIONS.keys():
                    passing_gene_count += 1 
            
            try:
                percentage_lims_genes_above = passing_gene_count / len(lims_genes)
            except:
                percentage_lims_genes_above = 0
            
            if test or (percentage_lims_genes_above >= MIN_PERCENT_LOCI_COVERED):
                if TNGS:
                    pncA_mutations = self.DF_LABORATORIAN[(self.DF_LABORATORIAN["tbprofiler_gene_name"] == "pncA")]
                    if "p.His57Asp" in pncA_mutations["tbprofiler_variant_substitution_aa"].tolist():
                        lineage.add("DNA of Mycobacterium bovis detected")
                    else:
                        lineage.add("DNA of Mycobacterium tuberculosis complex detected (not M. bovis)") 

                else:  
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
                lineage.add("DNA of Mycobacterium tuberculosis complex NOT detected")

            lineage = "; ".join(sorted(lineage))

        self.LINEAGE_ENGLISH = lineage
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

    def create_lims_report(self, TNGS: bool, MIN_PERCENT_LOCI_COVERED: float, OPERATOR: str) -> None:
        """Creates the LIMS report using the laboratorian pd.DataFrame
        
        Args:
            TNGS (bool): flag to indicate if the input data is tNGS or not
            MIN_PERCENT_LOCI_COVERED (float): the minimum percentage of the loci that need to be covered to call a lineage
            OPERATOR (str): the operator who ran the sequencing
        
        Returns:
            None
        """
        lims_report = {
            "Sample_Name": self.SAMPLE_NAME, 
            "Lineage_ID": self.get_id(TNGS, MIN_PERCENT_LOCI_COVERED)
        }

        qc_mask = ~self.DF_LABORATORIAN.apply(
            lambda row: row["tbprofiler_variant_substitution_nt"] in self.POSITIONAL_QC_FAILS.get(row["tbprofiler_gene_name"], set()), axis=1
        )
        qc_pass_laboratorian = self.DF_LABORATORIAN[qc_mask]

        self.logger.debug("LIMS:create_lims_report:Now iterating through each LIMS antimicrobial code")
        for drug, drug_gene_dictionary in globals_.DRUG_COLUMNS_TO_GENE_COLUMNS.items():
            for antimicrobial_code, gene_codes in drug_gene_dictionary.items():
                
                drug_associated_rows = qc_pass_laboratorian[qc_pass_laboratorian["antimicrobial"] == drug]
                gene_drug_associated_rows = drug_associated_rows.loc[drug_associated_rows["tbprofiler_gene_name"].isin(gene_codes.keys())]
                
                potential_mdl_resistances = gene_drug_associated_rows["mdl_interpretation"].tolist()
                try:
                    max_mdl_resistance = max(potential_mdl_resistances, key=lambda x: self.RESISTANCE_RANKING[x])
                except ValueError:
                    # no mutations associated with this drug -- this should only happen if all mutations were masked due to failing QC; set to "Insufficient Coverage"
                    max_mdl_resistance = "Insufficient Coverage"
                
                self.logger.debug("LIMS:create_lims_report:The max MDL resistance for this antimicrobial ({}) is {}".format(drug, max_mdl_resistance))
                    
                lims_report[antimicrobial_code] = self.convert_annotation(max_mdl_resistance, drug) 

                for gene, gene_code in gene_codes.items():
                    gene_subset = gene_drug_associated_rows[gene_drug_associated_rows["tbprofiler_gene_name"] == gene]
                    mutation_list = []
                    # the antimicrobial code message may need to change depending on the mutations present
                    # additionally, we must now provide content for the gene codes associated with this drug
                    if max_mdl_resistance in ["WT", "Insufficient Coverage", "NA"]:
                        # if the gene has insufficient coverage without a valid deletion (NOT in GENES_WITH_VALID_DELETIONS but in LOW_DEPTH_OF_COVERAGE_LIST), set to "No sequence"
                        if gene in self.LOW_DEPTH_OF_COVERAGE_LIST and gene not in self.GENES_WITH_VALID_DELETIONS.keys():
                            lims_report[gene_code] = "No sequence"
                            lims_report[antimicrobial_code] = "Pending Retest"
                        else:
                            lims_report[gene_code] = "No mutations detected"
                        
                    elif max_mdl_resistance in ["S"]:
                        if gene == "rpoB" and drug == "rifampicin":
                            # by default, assume no RRDR synonymous mutations
                            lims_report[antimicrobial_code] = "Predicted susceptibility to rifampicin"

                            # check if there is a mutation appears in the RRDR (codons 426-452) (in the special-positions dictionary)
                            # reduce gene_subset to only those that are synonymous type
                            synonymous_subset = gene_subset[gene_subset["tbprofiler_variant_substitution_type"] == "synonymous_variant"]
                               
                            for aa_mutation in synonymous_subset["tbprofiler_variant_substitution_aa"].tolist():
                                position_aa = globals_.get_position(aa_mutation)
                                
                                if globals_.is_mutation_within_range(position_aa, self.SPECIAL_POSITIONS[gene]):
                                    # synonymous RRDR mutation - update antimicrobial code message
                                    lims_report[antimicrobial_code] = "Predicted susceptibility to rifampicin. The detected synonymous mutation(s) do not confer resistance"
                                    mutation_list.append("{} [synonymous]".format(aa_mutation))
                            
                            lims_report[gene_code] = "; ".join(mutation_list) if len(mutation_list) > 0 else "No high confidence mutations detected"

                        else: 
                            if len(gene_subset) > 0 and max(gene_subset["mdl_interpretation"].tolist(), key=lambda x: self.RESISTANCE_RANKING[x]) == "S":
                                lims_report[gene_code] = "No high confidence mutations detected"
                            else:
                                # no mutations for this gene, wildtype
                                lims_report[gene_code] = "No mutations detected"
                    
                    elif max_mdl_resistance in ["R", "U"]:
                        # only report R or U mutations (and RRDR synonymous mutations that are S)
                        nt_list = gene_subset["tbprofiler_variant_substitution_nt"].tolist()
                        aa_list = gene_subset["tbprofiler_variant_substitution_aa"].tolist()
                        type_list = gene_subset["tbprofiler_variant_substitution_type"].tolist()
                        mdl_list = gene_subset["mdl_interpretation"].tolist()
        
                        for mutation, aa_mutation, mutation_type, mdl_interpretation in zip(nt_list, aa_list, type_list, mdl_list):
                            if ((mdl_interpretation in ["R", "U"]) 
                                or (gene == "rpoB" 
                                    and globals_.is_mutation_within_range(globals_.get_position(aa_mutation), self.SPECIAL_POSITIONS[gene]) 
                                    and mutation_type == "synonymous_variant" 
                                    and mdl_interpretation == "S")):
                                
                                # the following if only captures synonymous mutations if rpoB RRDR
                                if mutation_type == "synonymous_variant":
                                    substitution = "{} [synonymous]".format(substitution)
                                # R or U mutation
                                # report only amino acid mutations unless not applicable/blank/p.0?, in which case report nucleotide mutation
                                elif aa_mutation != "" and aa_mutation != "p.0?" and aa_mutation != "NA": 
                                    substitution = "{}".format(aa_mutation)
                                else:
                                    substitution = "{}".format(mutation)

                                mutation_list.append(substitution)
                            
                        if gene == "rpoB": # check special rpob positions to determine if we need to change the antimicrobial code message
                            if max_mdl_resistance == "R":
                                rpob_specific_mutations_counter = 0

                                for mutation in mutation_list:
                                    if any(rpob_mutation in mutation for rpob_mutation in self.RPOB_MUTATIONS):
                                        rpob_specific_mutations_counter += 1
                                    else:
                                        if mdl_interpretation == "R":
                                            # if there is an "R" mutation NOT in the special list, stop checking since we want the regular message here
                                            rpob_specific_mutations_counter = 0
                                            break

                                if rpob_specific_mutations_counter > 0:
                                    lims_report[antimicrobial_code] = "Predicted low-level resistance to rifampicin. May test susceptible by phenotypic methods."
                                else:
                                    lims_report[antimicrobial_code] = "Predicted resistance to rifampicin"


                    if len(mutation_list) > 0:
                        lims_report[gene_code] = "; ".join(mutation_list)
                    elif ((gene_code not in lims_report.keys()) 
                          or (lims_report[gene_code] not in ["No sequence", "No mutations detected"])):
                        # if max_mdl_resistance for this gene is "S" then set to "No high confidence mutations detected"
                        if len(gene_subset) > 0 and max(gene_subset["mdl_interpretation"].tolist(), key=lambda x: self.RESISTANCE_RANKING[x]) == "S":
                                lims_report[gene_code] = "No high confidence mutations detected"
                        else:
                            # no mutations for this gene, wildtype
                            lims_report[gene_code] = "No mutations detected"
                                        
        DF_LIMS = pd.DataFrame([lims_report])
        DF_LIMS["Analysis_Date"]= datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        DF_LIMS["Operator"] = OPERATOR
        DF_LIMS["Lineage"] = self.LINEAGE
        
        # write to file
        DF_LIMS.to_csv("{}.lims_report.csv".format(self.OUTPUT_PREFIX), index=False)
        DF_LIMS.T.to_csv("{}.lims_report.transposed.csv".format(self.OUTPUT_PREFIX), header=False)
        
        self.logger.info("LIMS:create_lims_report:LIMS report created, now exiting function\n")
