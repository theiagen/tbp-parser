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

  def __init__(self, logger, input_json, output_prefix, tngs):
    self.logger = logger
    self.input_json = input_json
    self.output_prefix = output_prefix
    self.tngs = tngs
  
  def get_id(self, test=False):
    """
    Returns the lineage in English for LIMS
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
      
      self.logger.debug("LIMS:Calculating the percentage of LIMS genes above the coverage threshold")
      percentage_limit = 0.7
      
      try:
        if self.tngs:
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
    
      if test or (percentage_lims_genes_above >= percentage_limit):
        if self.tngs:
          self.logger.debug("LIMS:[tNGS only] The sequencing method is tNGS; now checking for a His57Asp mutation in pncA")
          pncA_mutations = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == "pncA")]
          if "p.His57Asp" in pncA_mutations["tbprofiler_variant_substitution_aa"].tolist():
            self.logger.debug("LIMS:[tNGS only] p.His57Asp detected in pncA, lineage is likely M. bovis")
            lineage.add("DNA of Mycobacterium bovis detected")
          else:
            self.logger.debug("LIMS:[tNGS only] p.His57Asp not detected in pncA, lineage is likely M. tuberculosis")
            lineage.add("DNA of Mycobacterium tuberculosis complex detected (not M. bovis)") 
          
        else:  
          self.logger.debug("LIMS:The sequencing method is WGS AND the percentage of LIMS genes is GREATER than {}% ({}); now checking the TBProfiler lineage calls".format(percentage_limit * 100, percentage_lims_genes_above))

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
        self.logger.debug("LIMS:The percentage of LIMS genes is LESS than {}% ({}); assuming NOT M.tb".format(percentage_limit * 100, percentage_lims_genes_above))
        lineage.add("DNA of Mycobacterium tuberculosis complex NOT detected")
       
      lineage = "; ".join(sorted(lineage))
        
    globals_.LINEAGE_ENGLISH = lineage
    self.logger.debug("LIMS:The LIMS ID is: '{}'".format(globals_.LINEAGE_ENGLISH))
    self.logger.debug("LIMS:The TBProfiler lineage is: '{}'".format(globals_.LINEAGE))
    self.logger.info("LIMS:Finished getting lineage and ID, now exiting function")
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
   
  def apply_lims_rules(self, gene_dictionary, DF_LIMS, max_mdl_resistance, antimicrobial_code, responsible_gene):
    """
    This function implements several parsing rules for the LIMS report.
    Explanation of input variables:
    - gene_dictionary = the genes matched to their associated LIMS codes
    - DF_LIMS = the LIMS dataframe
    - antimicrobial_code = the LIMS code for the drug
    """ 
    self.logger.info("LIMS:Within LIMS class apply_lims_rules function; applying rules to mutations associated with {}".format(globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]))
    
    antimicrobial_name = globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]
    mutations_per_gene = {}
    
    # get all MDL interpretations for the responsible gene(s) for the max mdl interpretations in case of recalculation
    all_responsible_mdl_interpretations = {}
    for responsible_gene_name in responsible_gene:
      all_responsible_mdl_interpretations[responsible_gene_name] = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == responsible_gene_name) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["mdl_interpretation"].tolist()
      
    self.logger.debug("LIMS:The MDL interpretations for the responsible gene(s) are: {}".format(all_responsible_mdl_interpretations))  
      
    for gene, gene_code in gene_dictionary.items():
      DF_LIMS[gene_code] = ""
      non_s_mutations = 0
      
      # check to see if we are reporting on this gene for the LIMS report
      if gene in globals_.GENES_FOR_LIMS:
        # get all mutations, their types, interpretations, and warnings associated with a particular gene and it's associated antimicrobial
        nt_mutations_per_gene = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["tbprofiler_variant_substitution_nt"].tolist()
        aa_mutations_per_gene = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["tbprofiler_variant_substitution_aa"].tolist()
        mutation_types_per_gene = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["tbprofiler_variant_substitution_type"].tolist()
        mdl_interpretations = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["mdl_interpretation"].tolist()
        warnings = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["warning"].tolist()
        read_supports = globals_.DF_LABORATORIAN[(globals_.DF_LABORATORIAN["tbprofiler_gene_name"] == gene) & (globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial_name)]["read_support"].tolist()  
                 
        self.logger.debug("LIMS:The following mutations belong to this gene ({}) and are associated with this drug ({})".format(gene, antimicrobial_name))
        self.logger.debug("LIMS:Amino acid mutations: {}".format(aa_mutations_per_gene))
        self.logger.debug("LIMS:Nucleotide mutations: {}".format(nt_mutations_per_gene))
        self.logger.debug("LIMS:Their corresponding MDL interpretations: {}".format(mdl_interpretations))
        
        # check if there are any matching amino acid positions;
        # if so, we want to keep the row with the higher read support
        self.logger.debug("LIMS:Considering if any mutations have identical amino acid positions and keeping only the one with higher read support")
        removal_list = set()
        for mutation in aa_mutations_per_gene:
          current_index = aa_mutations_per_gene.index(mutation)
          
          # get a list of all other mutations for this gene except the current index for comparison          
          aa_positions_original = {aa_mutation:globals_.get_position(aa_mutation) for i, aa_mutation in enumerate(aa_mutations_per_gene) if i != current_index}
          aa_positions_flattened = {aa_mutation:subposition for aa_mutation, position in aa_positions_original.items() for subposition in position}
          
          if mutation not in ["Insufficient Coverage", "NA", "WT"]:
            other_positions = set(aa_positions_flattened.values())
            current_positions = set(globals_.get_position(mutation))
            same_positions = list(current_positions & other_positions)
            
            if len(same_positions) > 0:
              matching_index = aa_mutations_per_gene.index(list(aa_positions_flattened.keys())[list(aa_positions_flattened.values()).index(same_positions[0])])
              if matching_index == current_index:
                self.logger.debug("LIMS:We found the same mutation, checking for the second one instead")
                matching_index = aa_mutations_per_gene.index(list(aa_positions_flattened.keys())[list(aa_positions_flattened.values()).index(same_positions[0])], current_index + 1)

              # confirm the mutation is actually the same before carrying on
              if aa_mutations_per_gene[matching_index] == aa_mutations_per_gene[current_index]:
                self.logger.debug("LIMS:The current mutation has a matching amino acid position to an additional mutation, now testing read support")
                if read_supports[current_index] > read_supports[matching_index]:
                  self.logger.debug("LIMS:The current mutation ({}) has higher read support ({}) than the other mutation ({}; {})".format(mutation, read_supports[current_index], aa_mutations_per_gene[matching_index], read_supports[matching_index]))
                  removal_list.add(aa_mutations_per_gene[matching_index]) # avoid removing items while iterating through object
                else:
                  self.logger.debug("LIMS:The other mutation ({}) has higher read support ({}) than the current mutation ({}; {})".format(aa_mutations_per_gene[matching_index], read_supports[matching_index], mutation, read_supports[current_index]))
                  removal_list.add(mutation)
          
          if ("This mutation is outside the expected region" in warnings[current_index]):
            removal_list.add(mutation)
            
        # remove all mutations that have lower read support
        if len(removal_list) > 0:
          for mutation in removal_list:
            if mutation in aa_mutations_per_gene:
              removal_index = aa_mutations_per_gene.index(mutation)
              aa_mutations_per_gene.pop(removal_index)
              nt_mutations_per_gene.pop(removal_index)
              mutation_types_per_gene.pop(removal_index)
              mdl_interpretations.pop(removal_index)
              warnings.pop(removal_index)
              read_supports.pop(removal_index)  
        
        # format all mutations that are associated with the drug appropriately
        for mutation in nt_mutations_per_gene:
          
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
            if globals_.RESISTANCE_RANKING[mdl_interpretations[index]] >= 2: 
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
              
              if (max([globals_.RESISTANCE_RANKING[interpretation] for gene_set in all_responsible_mdl_interpretations.values() for interpretation in gene_set]) != globals_.RESISTANCE_RANKING[max_mdl_resistance[0]]) and gene in responsible_gene:
                max_mdl_resistance = [annotation for annotation, rank in globals_.RESISTANCE_RANKING.items() if rank == max([globals_.RESISTANCE_RANKING[interpretation] for gene_set in all_responsible_mdl_interpretations.values() for interpretation in gene_set])]
                self.logger.debug("LIMS:The maximum needed to be reevaluated; the potential new max_mdl_resistance is now {}".format(max_mdl_resistance[0]))
                
                if DF_LIMS[antimicrobial_code][0] != "Pending Retest":
                  self.logger.debug("LIMS:Now changing the antimicrobial code value for this gene since max_mdl_resistance changed")
                  DF_LIMS[antimicrobial_code] = self.convert_annotation(max_mdl_resistance[0], antimicrobial_name)   
                else:
                  self.logger.debug("LIMS:The antimicrobial code value is currently 'Pending Retest'; this shouldn't be overwritten.")
              else:
                self.logger.debug("LIMS:The maximum did not need to be reevaluated; the max_mdl_resistance is still {}".format(max_mdl_resistance[0]))
          # the mutation is of decent quality and non-S, we want to report all non-synonymous mutations UNLESS rpoB RRDR (see Variant l.145 for explanation)
          elif ((mutation_type != "synonymous_variant" and mdl_interpretations[index] != "S") 
                or (gene == "rpoB" and globals_.is_within_range(position_aa, globals_.SPECIAL_POSITIONS[gene]))):
              if aa_mutation != "" and aa_mutation != "p.0?": # report only amino acid mutations unless not applicable/blank, in which case report nucleotide mutation
                # and also p.0? and p.0
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
                if any(rpob_mutation in mutation for rpob_mutation in globals_.RPOB_MUTATIONS):
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
          maximum_ranking = max([globals_.RESISTANCE_RANKING[interpretation] for interpretation in mdl_interpretations])
          # see the globals_.RESISTANCE_RANKING dictionary for the ranking of each MDL interpretation
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
          self.logger.debug("LIMS:There are no mutations for this gene ({}) associated with this drug ({})".format(gene, antimicrobial_name))
          DF_LIMS[gene_code] = "No mutations detected"

        if (DF_LIMS[gene_code][0] == "No sequence" or "Insufficient Coverage" in mdl_interpretations) and max_mdl_resistance[0] != "R" and gene not in globals_.GENES_WITH_DELETIONS:
          self.logger.debug("LIMS:This gene ({}) has insufficient coverage without a valid deletion and no other mutations associated with this antimicrobial ({}) are 'R'; changing antimicrobial output to 'Pending Retest'".format(gene, antimicrobial_name))
          DF_LIMS[antimicrobial_code] = "Pending Retest"    

    self.logger.info("LIMS:Finished applying rules to mutations associated with {}, exiting function".format(globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]))
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
    self.logger.info("LIMS:Within LIMS class create_lims_report function")
    DF_LIMS = pd.DataFrame({
      "MDL sample accession numbers": globals_.SAMPLE_NAME, 
      "M_DST_A01_ID": self.get_id()
      }, index=[0])
    
    self.logger.debug("LIMS:Now iterating through each LIMS antimicrobial code")
    for antimicrobial_code, gene_dictionary in globals_.ANTIMICROBIAL_CODE_TO_GENES.items():
      drug_name = globals_.ANTIMICROBIAL_CODE_TO_DRUG_NAME[antimicrobial_code]
      # get the MDL interpretations for all genes **FOR THE LIMS REPORT** associated with this drug             
      potential_mdl_resistances = globals_.DF_LABORATORIAN[globals_.DF_LABORATORIAN["antimicrobial"] == drug_name].loc[globals_.DF_LABORATORIAN["tbprofiler_gene_name"].isin(gene_dictionary.keys())]
      
      # remove any interpretations with failed quality warnings
      potential_mdl_resistances = potential_mdl_resistances.loc[~potential_mdl_resistances["warning"].str.contains("Failed quality in the mutation position")]
      potential_mdl_resistances = potential_mdl_resistances["mdl_interpretation"].tolist()

      # initalize list of genes responsible for the max resistance
      responsible_gene = set()
      
      # get the maximum resistance for the drug
      try:
        max_mdl_resistance = [annotation for annotation, rank in globals_.RESISTANCE_RANKING.items() if rank == max([globals_.RESISTANCE_RANKING[interpretation] for interpretation in potential_mdl_resistances])]
        responsible_gene = set(globals_.DF_LABORATORIAN[globals_.DF_LABORATORIAN["antimicrobial"] == drug_name].loc[globals_.DF_LABORATORIAN["mdl_interpretation"] == max_mdl_resistance[0]]["tbprofiler_gene_name"].tolist())

        # keep only LIMS genes
        responsible_gene = responsible_gene.intersection(set(gene_dictionary.keys()))
        self.logger.debug("LIMS:The gene(s) responsible for the max MDL resistance for this antimicrobial ({}) is/are {}".format(drug_name, responsible_gene))
      except:
        max_mdl_resistance = ["NA"]
      self.logger.debug("LIMS:The max MDL resistance for this antimicrobial ({}) is {}".format(drug_name, max_mdl_resistance[0]))

      DF_LIMS[antimicrobial_code] = self.convert_annotation(max_mdl_resistance[0], drug_name) 

      DF_LIMS = self.apply_lims_rules(gene_dictionary, DF_LIMS, max_mdl_resistance, antimicrobial_code, responsible_gene)
    
    DF_LIMS["Analysis date"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    DF_LIMS["Operator"] = globals_.OPERATOR
    
    # add lineage
    DF_LIMS["M_DST_O01_Lineage"] = globals_.LINEAGE

    # write to file
    DF_LIMS.to_csv("{}.lims_report.csv".format(self.output_prefix), index=False)
    self.logger.info("LIMS:LIMS report created, now exiting function\n")