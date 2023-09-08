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
   
  def apply_lims_rules(self, gene_dictionary, max_mdl_resistance, mutations_per_gene, DF_LIMS, non_s_mutations, antimicrobial_code):
    """
    This function implements several parsing rules for the LIMS report.
    """
    for gene, gene_code in gene_dictionary.items():

      if gene in globals.GENES_FOR_LIMS:
        nt_mutations_per_gene = globals.DF_LABORATORIAN[globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]["tbprofiler_variant_substitution_nt"].unique().tolist()
        aa_mutations_per_gene = globals.DF_LABORATORIAN[globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]["tbprofiler_variant_substitution_aa"].unique().tolist()
        mutation_types_per_gene = globals.DF_LABORATORIAN[globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]["tbprofiler_variant_substitution_type"].unique().tolist()
        mdl_interpretation = globals.DF_LABORATORIAN[globals.DF_LABORATORIAN["tbprofiler_gene_name"] == gene]["mdl_interpretation"].tolist()
        
        if max_mdl_resistance[0] == "S" and gene != "rpoB":
          mutations_per_gene[gene] = "No high confidence mutations detected"
        else:         
          # format all mutations associated with a particular antimicrobial appropriately
          for mutation in nt_mutations_per_gene:
            index = nt_mutations_per_gene.index(mutation)
            
            # perform some data clean-up:
            if mutation == "WT":
              mutation = ""
            
            try:
              aa_mutation = aa_mutations_per_gene[index]
            except:
              aa_mutation = ""
              
            if aa_mutation == "NA" or aa_mutation == "WT":
              aa_mutation = "" 
            if gene == "rpoB":           
              aa_mutation_position = globals.get_position(aa_mutation)
              
            try:  
              mutation_type = mutation_types_per_gene[index]
            except: 
              mutation_type = mutation_types_per_gene[0]
            
            # report all non-synonymous mutations UNLESS rpoB RRDR
            if mutation_type != "synonymous_variant" or (gene == "rpoB" and (globals.SPECIAL_POSITIONS[gene][1] <= aa_mutation_position <= globals.SPECIAL_POSITIONS[gene][2])):
              substitution = "{} ({})".format(mutation, aa_mutation)
              
              # the following if only captures synonymous mutations if rpoB RRDR mutations
              if mutation_type == "synonymous_variant":
                substitution = "{} [synonymous]".format(mutation)
              
              # add to dictionary
              if gene not in mutations_per_gene.keys():
                mutations_per_gene[gene] = [substitution]
              else:
                mutations_per_gene[gene] = "{}; {}".format("".join(mutations_per_gene[gene]), substitution)
            else:
              mutations_per_gene[gene] = "No high confidence mutations detected"
              
          # if that gene has mutations associated with it, perform additional filtration -- also check for coverage???
          if gene in mutations_per_gene.keys():
            DF_LIMS[gene_code] = mutations_per_gene[gene]
            if mdl_interpretation[index] != "S" or mdl_interpretation[index] != "U" or mdl_interpretation[index] != "WT":
              non_s_mutations += 1

            if gene == "rpoB":
              if mdl_interpretation[index] == "R":
                non_rpob_specific_mutations_counter = 0
                
                for mutation in mutations_per_gene[gene]:
                  if ~ any(rpob_mutation in mutation for rpob_mutation in globals.RPOB_MUTATIONS):
                    non_rpob_specific_mutations_counter += 1
                
                if non_rpob_specific_mutations_counter == 0:
                  DF_LIMS[antimicrobial_code] = "Predicted low-level resistance to rifampin. May test susceptible by phenotypic methods."
                else:
                  DF_LIMS[antimicrobial_code] = "Predicted resistance to rifampin"

              if mdl_interpretation[index] == "S":
                if len(mutations_per_gene[gene]) > 0: 
                  non_synonymous_count = 0
                  
                  for mutation in mutations_per_gene[gene]: 
                    # if any NON synonymous mutations were identified, we want to keep the original output
                    if "synonymous" not in mutation: 
                      non_synonymous_count += 1 
                  # otherwise, the only synonymous mutations were identified in rpoB RRDR
                  if non_synonymous_count == 0: 
                    DF_LIMS[antimicrobial_code] = "No mutations associated with resistance to rifampin detected. The detected synonymous mutation(s) do not confer resistance but may result in false-resistance in PCR-based assays targeting the rpoB RRDR."

            elif mdl_interpretation[index] == "S" and non_s_mutations == 0:
              DF_LIMS[gene_code] = "No high confidence mutations detected"
            elif mdl_interpretation[index] == "WT":
              DF_LIMS[gene_code] = "No mutations detected"
            elif mdl_interpretation[index] == "Insufficient Coverage":
              DF_LIMS[gene_code] = "No sequence"
        
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

      DF_LIMS[antimicrobial_code] = self.convert_annotation(max_mdl_resistance[0], drug_name)      

      mutations_per_gene = {}
      non_s_mutations = 0

      DF_LIMS = self.apply_lims_rules(gene_dictionary, max_mdl_resistance, mutations_per_gene, DF_LIMS, non_s_mutations, antimicrobial_code)
           
    DF_LIMS["Analysis date"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    DF_LIMS["Operator"] = globals.OPERATOR
    
    # write to file
    DF_LIMS.to_csv("{}.lims_report.csv".format(self.output_prefix), index=False)