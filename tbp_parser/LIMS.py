import globals
import pandas as pd
import json
import datetime

class LIMS:
  """
  This class creates the CDPH LIMS report.
  
  It has four functions:
    - get_lineage: returns the lineage in English for LIMS
    - convert_annotation: converts the resistance annotation and the target drug
      into the LIMS language
    - get_mutation_position: returns the position where a mutation occurs
    - apply_lims_rules: implements several parsing rules for the LIMS report
    - create_lims_report: creates the LIMS report CSV file
  """

  def __init__(self, logger, input_json, output_prefix):
    self.logger = logger
    self.input_json = input_json
    self.output_prefix = output_prefix
  
  def get_lineage(self):
    """
    Returns the lineage in English for LIMS
    """
    self.logger.info("Within LIMS class get_lineage function")
    
    self.logger.debug("Calculating the percentage of genes about the coverage threshold")
    number_of_genes_above_coverage_threshold = sum(int(value) >= globals.COVERAGE_THRESHOLD for value in globals.COVERAGE_DICTIONARY.values())
    percentage_above = number_of_genes_above_coverage_threshold / len(globals.COVERAGE_DICTIONARY)
    self.logger.debug("The percentage of genes above the coverage threshold is {}".format(percentage_above))
    
    # if the percentage of genes above the coverage threshold is greater than 90%, then we can call the lineage
    percentage_limit = 0.9
    
    with open(self.input_json) as json_fh:
      input_json = json.load(json_fh)
      
      # set default value for lineage
      lineage = "DNA of Mycobacterium tuberculosis complex detected (not M. bovis and not M. tb)"
      detected_lineage = input_json["main_lin"]
      
      if detected_lineage == "":
        if percentage_above >= percentage_limit:
          lineage = "DNA of Mycobacterium tuberculosis complex detected"
        else:
           lineage = "DNA of Mycobacterium tuberculosis complex NOT detected"
      elif "lineage" in detected_lineage:
        lineage = "DNA of Mycobacterium tuberculosis species detected"
      elif "BCG" in detected_lineage:
        lineage = "DNA of Mycobacterium bovis BCG detected"
      elif "bovis" in detected_lineage:
        lineage = "DNA of non-BCG Mycobacterium bovis detected"
      elif "La1" in detected_lineage:
        lineage = "DNA of M. tuberculosis complex detected (M. bovis)"
      
    self.logger.debug("The lineage is: {}".format(lineage))
    self.logger.info("Finished getting lineage, now exiting function")
    return lineage
  
  def convert_annotation(self, annotation, drug):
    """
    This function takes the resistance annotation and the target drug
    and converts it into the LIMS language
    """
    message = "No mutations associated with resistance to {} detected".format(drug)
    if annotation == "R":
      message = "Mutation(s) associated with resistance to {} detected".format(drug)
    elif (annotation == "R-Interim") or (annotation == "U"):
      message = "The detected mutation(s) have uncertain significance. Resistance to {} cannot be ruled out".format(drug)
    elif annotation == "Insufficient Coverage":
      message = "Pending Retest"      
      
    return message
   
  def apply_lims_rules(self, gene_dictionary, DF_LIMS, antimicrobial_code):
    """
    This function implements several parsing rules for the LIMS report.
    Explanation of input variables:
    - gene_dictionary = the genes matched to their associated LIMS codes
    - DF_LIMS = the LIMS dataframe
    - antimicrobial_code = the LIMS code for the drug
    """ 
    self.logger.info("Within LIMS class apply_lims_rules function; applying rules to mutations associated with {}".format(globals.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]))
    
    antimicrobial_name = globals.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]
    mutations_per_gene = {}

    for gene, gene_code in gene_dictionary.items():
      non_s_mutations = 0
      # check to see if we are reporting on this gene for the LIMS report
      if gene in globals.GENES_FOR_LIMS:
        # get all mutations, their types, interpretations, and warnings associated with a particular gene and it's associated antimicrobial
        nt_mutations_per_gene = globals.DF_LABORATORIAN[(globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["tbprofiler_variant_substitution_nt"].tolist()
        aa_mutations_per_gene = globals.DF_LABORATORIAN[(globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["tbprofiler_variant_substitution_aa"].tolist()
        mutation_types_per_gene = globals.DF_LABORATORIAN[(globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["tbprofiler_variant_substitution_type"].tolist()
        mdl_interpretations = globals.DF_LABORATORIAN[(globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["mdl_interpretation"].tolist()
        warnings = globals.DF_LABORATORIAN[(globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["warning"].tolist()
        
        self.logger.debug("The following mutations belong to this gene ({}) and are associated with this drug ({})".format(gene, antimicrobial_name))
        self.logger.debug("Nucleotide mutations: {}".format(nt_mutations_per_gene))
        self.logger.debug("Their corresponding MDL interpretations: {}".format(mdl_interpretations))
        
        # format all mutations that are associated with the drug appropriately
        for mutation in nt_mutations_per_gene:
          
          # grab the index so we can also grab the other information associated with this mutation
          index = nt_mutations_per_gene.index(mutation)
          aa_mutation = aa_mutations_per_gene[index]
          mutation_type = mutation_types_per_gene[index]
          
          # if it is rpoB, we want to get the position of the mutation to check if it is in RRDR
          if gene == "rpoB":           
            aa_mutation_position = globals.get_position(aa_mutation)
          
          # convert WT and Insufficient Coverage to empty strings
          if mutation == "WT" or mutation == "Insufficient Coverage" or mutation == "NA":
            mutation = ""
            
          # convert NA,  WT, and Insufficient Coverage to empty strings
          if aa_mutation == "NA" or aa_mutation == "WT" or aa_mutation == "Insufficient Coverage":
            aa_mutation = ""
              
          # do not add the mutation if the particular mutation has low quality or is blank
          if ("Failed quality in the mutation position" in warnings[index]) or (mutation == ""):
            self.logger.debug("This mutation (\"{}\", origin gene: {}) is not being added to the LIMS report because it failed quality in the mutation position, was WT, or had insufficient locus coverage".format(mutation, gene))
            
          # the mutation is of decent quality and non-S, we want to report all non-synonymous mutations UNLESS rpoB RRDR
          elif (mutation_type != "synonymous_variant" and mdl_interpretations[index] != "S") or (gene == "rpoB" and (globals.SPECIAL_POSITIONS[gene][0] <= aa_mutation_position <= globals.SPECIAL_POSITIONS[gene][1])):
            substitution = "{} ({})".format(mutation, aa_mutation)
            
            # the following if only captures synonymous mutations if rpoB RRDR mutations
            if mutation_type == "synonymous_variant":
              substitution = "{} [synonymous]".format(substitution)
            
            # add the formatted mutation to mutations_per_gene dictionary (formatted as a list)
            if gene not in mutations_per_gene.keys():
              mutations_per_gene[gene] = [substitution]
            elif substitution not in mutations_per_gene[gene]:
                mutations_per_gene[gene] = "{}; {}".format("".join(mutations_per_gene[gene]), substitution)
          else:
            self.logger.debug("This mutation (\"{}\", origin gene: {}) is not being added to the LIMS report because it is a non-rpoB RRDR \"S\" mutation".format(mutation, gene))
        
        # Mutations for a particular gene have been added to the mutations_per_gene dictionary.
        # if that gene has mutations associated with it, we want to perform some additional filtration,
        # especially if it is RRDR rpoB
        if gene in mutations_per_gene.keys():
          self.logger.debug("Adding the mutations ({}) associated with this gene ({}) to the DF_LIMS dataframe".format(mutations_per_gene[gene], gene))
          DF_LIMS[gene_code] = mutations_per_gene[gene]
        
          if gene == "rpoB":
            if mdl_interpretations[index] == "R":
              self.logger.debug("This gene is rpoB, now checking to see if all of the mutations belong to the special position list")
              rpob_specific_mutations_counter = 0
              
              for mutation in mutations_per_gene[gene]:
                if any(rpob_mutation in mutation for rpob_mutation in globals.RPOB_MUTATIONS):
                  rpob_specific_mutations_counter += 1
              
              if rpob_specific_mutations_counter > 0:
                self.logger.debug("The only rpoB mutations are in the special rpoB mutation list, changing the output message")
                DF_LIMS[antimicrobial_code] = "Predicted low-level resistance to rifampin. May test susceptible by phenotypic methods."
              else:
                self.logger.debug("Not all mutations are in the special rpoB mutation list, changing the output message")
                DF_LIMS[antimicrobial_code] = "Predicted resistance to rifampin"

            # if the maximum MDL resistance is S in this case, the mutation will always be an rpoB RRDR region one
            if mdl_interpretations[index] == "S":
              self.logger.debug("A synonymous mutations was identified in RRDR; changing output message")
              DF_LIMS[antimicrobial_code] = "Predicted susceptibility to rifampin. The detected synonymous mutation(s) do not confer resistance"

        elif gene == "rpoB":
          self.logger.debug("No mutations were identified in rpoB; changing output message")
          DF_LIMS[antimicrobial_code] = "Predicted susceptibility to rifampin"
            
        # make sure that there is a mutation associated with this gene-drug combo
        if len(nt_mutations_per_gene) > 0:
          maximum_ranking = max([globals.RESISTANCE_RANKING[interpretation] for interpretation in mdl_interpretations])
          # see the globals.RESISTANCE_RANKING dictionary for the ranking of each MDL interpretation
          # count up the number of non-S mutations; non-S mutations have a ranking > 2
          if maximum_ranking > 2:
            non_s_mutations += 1

          # change the gene_code to be something different depending on the MDL interpretations and/or number of non-s mutations
          if maximum_ranking == 2 and non_s_mutations == 0: # S
            DF_LIMS[gene_code] = "No high confidence mutations detected"
          elif maximum_ranking == 1: # WT
            DF_LIMS[gene_code] = "No mutations detected"
          elif maximum_ranking == 0 and "del" not in DF_LIMS[gene_code]: # Insufficient Coverage (and not a deletion)
            DF_LIMS[gene_code] = "No sequence"
        else:
          self.logger.debug("There are no mutations for this gene ({}) associated with this drug ({})".format(gene, antimicrobial_name))
          DF_LIMS[gene_code] = "No mutations detected"
    
    self.logger.info("Finished applying rules to mutations associated with {}, exiting function".format(globals.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]))
    return DF_LIMS
  
  def create_lims_report(self):
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
    self.logger.info("Within LIMS class create_lims_report function")
    DF_LIMS = pd.DataFrame({
      "MDL sample accession numbers": globals.SAMPLE_NAME, 
      "M_DST_A01_ID": self.get_lineage()
      }, index=[0])
    
    self.logger.debug("Now iterating through each LIMS antimicrobial code")
    for antimicrobial_code, gene_dictionary in globals.ANTIMICROBIAL_CODE_TO_GENES.items():
      drug_name = globals.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]
      
      potential_mdl_resistances = globals.DF_LABORATORIAN[globals.DF_LABORATORIAN["antimicrobial"] == drug_name]["mdl_interpretation"]
      
      # get the maximum resistance for the drug
      max_mdl_resistance = [annotation for annotation, rank in globals.RESISTANCE_RANKING.items() if rank == max([globals.RESISTANCE_RANKING[interpretation] for interpretation in potential_mdl_resistances])]
      self.logger.debug("The max MDL resistance for this antimicrobial ({}) is {}".format(drug_name, max_mdl_resistance[0]))

      DF_LIMS[antimicrobial_code] = self.convert_annotation(max_mdl_resistance[0], drug_name)            

      DF_LIMS = self.apply_lims_rules(gene_dictionary, DF_LIMS, antimicrobial_code)
           
    DF_LIMS["Analysis date"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    DF_LIMS["Operator"] = globals.OPERATOR
    
    # write to file
    DF_LIMS.to_csv("{}.lims_report.csv".format(self.output_prefix), index=False)
    self.logger.info("LIMS report created, now exiting function\n")