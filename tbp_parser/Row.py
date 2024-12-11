import globals

class Row() :
  """
  This class represents a row in the CDPH Laboratorian report.
  The __init__ function assigns each variant attribute to the appropriate
  column in the CDPH Laboratorian report; it also applies any coverage
  warnings if necessary.
  
  This class has five additional functions:
    - print: prints the row in a readable format
    - complete_row: finishes each row with the rest of the values needed
    - rank_annotation: ranks the WHO annotation based on resistance
    - annotation_to_LIMS: converts the WHO annotation and the target drug
      into returns the LIMS' report file appropriate annotation
    - describe_rationale: removes the 'noexpert' suffix in the case where
      the interpretation logic applied is not considered an expert rule
  """
  
  def __init__(self, logger, variant, who_confidence, drug, gene_name=None, depth=0, frequency=None, source="", tbdb_comment=""):
    self.logger = logger
    
    self.variant = variant
    self.who_confidence = who_confidence
    self.antimicrobial = drug

    # create a row for the sample if the gene is in the coverage dictionary (should occur for all except in tNGS where only certain genes are sequenced)
    if gene_name != "test":
      self.antimicrobial = self.antimicrobial.replace("rifampicin", "rifampin")

      self.sample_id = globals.SAMPLE_NAME

      # Initalizing the rest of the columns for the CDPH Laboratorian report
      # for when the variant is in the JSON file
      if variant is not None:
        self.logger.debug("ROW:Initalizing the Row object, the variant has been supplied.")
        
        try:
          self.tbprofiler_gene_name = self.variant.gene_name
        except:
          self.tbprofiler_gene_name = gene_name
          self.variant.gene_name = gene_name
        if self.tbprofiler_gene_name == "mmpR5":
          self.tbprofiler_gene_name = "Rv0678"
          
        self.logger.debug("ROW:The variant's gene name is {}".format(self.tbprofiler_gene_name))
        try:
          self.tbprofiler_locus_tag = self.variant.locus_tag
        except:
          self.tbprofiler_locus_tag = globals.GENE_TO_LOCUS_TAG[self.tbprofiler_gene_name]
        self.tbprofiler_variant_substitution_type = self.variant.type
        self.tbprofiler_variant_substitution_nt = self.variant.nucleotide_change
        self.tbprofiler_variant_substitution_aa = self.variant.protein_change
        self.logger.debug("ROW:This mutation is a {} with nucleotide change \"{}\" and protein change \"{}\"".format(self.tbprofiler_variant_substitution_type, self.tbprofiler_variant_substitution_nt, self.tbprofiler_variant_substitution_aa))
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
        
        # the following if statement only applies to rpoB
        # if self.tbprofiler_gene_name in globals.TNGS_REGIONS.keys():
        #   self.logger.debug("ROW:[tNGS only] This mutation's genomic position is outside the expected region.")
        #   self.warning.append("This mutation is outside the expected region")
        if self.tbprofiler_gene_name in globals.COVERAGE_DICTIONARY:
          if (self.depth < globals.MIN_DEPTH) or (float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) < globals.COVERAGE_THRESHOLD):
            self.logger.debug("ROW:The depth of coverage for this variant is {} and the coverage for the gene is {}; applying a locus warning".format(self.depth, globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]))
            if (float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) < globals.COVERAGE_THRESHOLD):
              globals.LOW_DEPTH_OF_COVERAGE_LIST.append(self.tbprofiler_gene_name)
              
              if "del" in self.tbprofiler_variant_substitution_nt or self.tbprofiler_gene_name in globals.GENES_WITH_DELETIONS:
                self.logger.debug("ROW:This is a deletion, no warning added for the locus unless it fails positional qc (checked next)")
              else:
                self.warning.append("Insufficient coverage in locus")
        else:
          self.logger.debug("ROW:This gene does not appear in the coverage dictionary. An additional warning will be given.")
          if self.tbprofiler_gene_name in globals.TNGS_REGIONS.keys():
            self.logger.debug("ROW:[tNGS only] This mutation's genomic position is outside the expected region")
            self.warning.append("This mutation is outside the expected region")
            self.logger.debug("ROW:[tNGS only] Rewriting this variant's interpretation to NA since it shouldn't exist")
            self.looker_interpretation = "NA"
            self.mdl_interpretation = "NA"
          else:
            self.warning.append("This mutation is outside the expected region")
        
        protein_position = globals.get_position(self.tbprofiler_variant_substitution_aa)
          
        # check to see if we need to apply a mutation warning 
        # (check rrs & rrl for low frequency and read support; 
        #  also check ethA & rpoB for specific protein position frequency)
        if ((self.depth < globals.MIN_DEPTH) 
             or ((self.tbprofiler_gene_name not in ["rrs", "rrl"] and 
                 (float(self.frequency) < globals.MIN_FREQUENCY or self.read_support < globals.MIN_READ_SUPPORT))
             or (self.tbprofiler_gene_name == "rrs" and 
                 (float(self.frequency) < globals.RRS_FREQUENCY or self.read_support < globals.RRS_READ_SUPPORT)) 
             or (self.tbprofiler_gene_name == "rrl" and 
                 (float(self.frequency) < globals.RRL_FREQUENCY or self.read_support < globals.RRL_READ_SUPPORT)) 
             or (self.tbprofiler_gene_name == "ethA" and 
                 237 in protein_position and float(self.frequency) < globals.ETHA237_FREQUENCY)
             or (self.tbprofiler_gene_name == "rpoB" and 
                 449 in protein_position and float(self.frequency) < globals.RPOB449_FREQUENCY))):
               
          if ("del" in self.tbprofiler_variant_substitution_nt and self.depth == 0 and self.read_support == 0):
            # placeholder
            self.logger.debug("ROW:This is an okay scenario")
          else:
            self.logger.debug("ROW:The depth of coverage for this variant is {}, the frequency is {}, and the read support is {}; applying an additional mutation position warning".format(self.depth, self.frequency, self.read_support))
            
            if self.tbprofiler_gene_name in globals.COVERAGE_DICTIONARY.keys():
              if ((float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) < globals.COVERAGE_THRESHOLD) and 
                  ("del" in self.tbprofiler_variant_substitution_nt or self.tbprofiler_gene_name in globals.GENES_WITH_DELETIONS)):
                self.logger.debug("ROW:This deletion failed in the mutation position and there was insufficient coverage locus, adding insufficient coverage warning")
                self.warning.append("Insufficient coverage in locus")
                  
            globals.MUTATION_FAIL_LIST.append(self.tbprofiler_variant_substitution_nt)
            self.warning.append("Failed quality in the mutation position")
          
        else:
          self.logger.debug("ROW:The depth of coverage for this variant is {}, the frequency is {}, and the read support is {}; no additional warning added for the mutation position".format(self.depth, self.frequency, self.read_support))
          if len(self.warning) == 0:
            self.warning = [""]        
          
        self.logger.debug("ROW:This variant has the following warnings: {}".format(self.warning))
      
        if "Failed quality in the mutation position" not in self.warning and "del" in self.tbprofiler_variant_substitution_nt:
          globals.GENES_WITH_DELETIONS.add(self.tbprofiler_gene_name)
          self.logger.debug("ROW:This is a deletion that passed positional qc, adding to set")

      # otherwise, the variant does not appear in the JSON file (or was not sequenced [tNGS]) and default NA/WT values need to be supplied
      else:
        self.logger.debug("ROW:Initializing the Row object, the variant has no information supplied. Defaulting to NA or WT values.")
       
        if gene_name == "mmpR5":
          self.tbprofiler_gene_name = "Rv0678"
        else: 
          self.tbprofiler_gene_name = gene_name
        try:
          self.tbprofiler_locus_tag = globals.GENE_TO_LOCUS_TAG[self.tbprofiler_gene_name]
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
        if (self.tbprofiler_gene_name in globals.COVERAGE_DICTIONARY.keys() 
            or self.tbprofiler_gene_name in globals.TNGS_REGIONS.keys()):
          try: 
            if float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) >= globals.COVERAGE_THRESHOLD:
              self.logger.debug("ROW:A warning does not need to be applied; setting the variant information to WT")
              self.tbprofiler_variant_substitution_type = "WT"
              self.tbprofiler_variant_substitution_nt = "WT"
              self.tbprofiler_variant_substitution_aa = "WT"
              self.looker_interpretation = "S"
              self.mdl_interpretation = "WT"
            else:
              self.logger.debug("ROW:A warning needs to be applied; setting the variant information to insufficient coverage")
              self.tbprofiler_variant_substitution_type = "Insufficient Coverage"
              self.tbprofiler_variant_substitution_nt = "WT"
              self.tbprofiler_variant_substitution_aa = "WT"
              self.looker_interpretation = "Insufficient Coverage"
              self.mdl_interpretation = "Insufficient Coverage"
              self.warning.append("Insufficient coverage in locus")
          except:
            self.logger.debug("ROW:[tNGS only]: The gene does not appear in the coverage dictionary, but is in the TNGS regions dictionary")
            self.logger.debug("ROW:[tNGS only]: This indicates that the gene was sequenced, but coverage was calculated under a different name")
            self.logger.debug("ROW:[tNGS only]: A coverage warning will be applied if at least one segment is under the coverage threshold")
            
            for segment in globals.TNGS_REGIONS[self.tbprofiler_gene_name]:
              if float(globals.COVERAGE_DICTIONARY[segment]) >= globals.COVERAGE_THRESHOLD:
                self.logger.debug("ROW:[tNGS only]: This segment has good coverage, checking the other segments.")
                self.tbprofiler_variant_substitution_type = "WT"
                self.tbprofiler_variant_substitution_nt = "WT"
                self.tbprofiler_variant_substitution_aa = "WT"
                self.looker_interpretation = "S"
                self.mdl_interpretation = "WT"
              else:
                self.logger.debug("ROW:[tNGS only]: This segment has poor coverage and a warning needs to be applied; setting the variant information to insufficient coverage")
                self.tbprofiler_variant_substitution_type = "Insufficient Coverage"
                self.tbprofiler_variant_substitution_nt = "WT"
                self.tbprofiler_variant_substitution_aa = "WT"
                self.looker_interpretation = "Insufficient Coverage"
                self.mdl_interpretation = "Insufficient Coverage"
                self.warning.append("Insufficient coverage in locus")
                break

        else:
          self.logger.debug("ROW:This gene ({}) was not sequenced".format(self.tbprofiler_gene_name))
          self.tbprofiler_variant_substitution_type = "NA"
          self.tbprofiler_variant_substitution_nt = "NA"
          self.tbprofiler_variant_substitution_aa = "NA"
          self.looker_interpretation = "NA"
          self.mdl_interpretation = "NA"
      
      try:
        self.gene_tier = globals.GENE_TO_TIER[self.tbprofiler_gene_name]
      except:
        try:
          # iterate through the tNGS regions to see if we have a match
          parent_genes = globals.TNGS_REGIONS.keys()
          parent_gene = [gene for gene in parent_genes if self.tbprofiler_gene_name in globals.TNGS_REGIONS[gene].keys()][0]
          self.logger.debug("ROW:[tNGS only]: The parent gene ({}) of this segment ({}) was identified; now adding tier".format(parent_gene, self.tbprofiler_gene_name))
          self.gene_tier = globals.GENE_TO_TIER[parent_gene]
          # now that coverage has been calculated, we can now rename the gene to be the parent gene name
          self.logger.debug("ROW:[tNGS only]: Renaming the gene segment ({}) to be the parent gene name ({})".format(self.tbprofiler_gene_name, parent_gene))
          self.tbprofiler_gene_name = parent_gene
    
        except:
          self.gene_tier = "NA"
      
      self.source = source
      self.tbdb_comment = tbdb_comment
                    
  def print(self):
    """
    This function prints the row in a readable format.
    """
    self.logger.debug("ROW:Now printing the row in a readable format:")
    self.logger.debug("ROW:\tsample_id: {}".format(self.sample_id))
    self.logger.debug("ROW:\ttbprofiler_gene_name: {}".format(self.tbprofiler_gene_name))
    self.logger.debug("ROW:\ttbprofiler_locus_tag: {}".format(self.tbprofiler_locus_tag))
    self.logger.debug("ROW:\ttbprofiler_variant_substitution_type: {}".format(self.tbprofiler_variant_substitution_type))
    self.logger.debug("ROW:\ttbprofiler_variant_substitution_nt: {}".format(self.tbprofiler_variant_substitution_nt))
    self.logger.debug("ROW:\ttbprofiler_variant_substitution_aa: {}".format(self.tbprofiler_variant_substitution_aa))
    self.logger.debug("ROW:\tconfidence: {}".format(self.confidence))
    self.logger.debug("ROW:\tantimicrobial: {}".format(self.antimicrobial))
    self.logger.debug("ROW:\tlooker_interpretation: {}".format(self.looker_interpretation))
    self.logger.debug("ROW:\tmdl_interpretation: {}".format(self.mdl_interpretation))
    self.logger.debug("ROW:\tdepth: {}".format(self.depth))
    self.logger.debug("ROW:\tfrequency: {}".format(self.frequency))
    self.logger.debug("ROW:\tread_support: {}".format(self.read_support))
    self.logger.debug("ROW:\trationale: {}".format(self.rationale))
    self.logger.debug("ROW:\twarning: {}".format(self.warning))
    self.logger.debug("ROW:\ttier: {}".format(self.gene_tier))
    self.logger.debug("ROW:\tsource: {}".format(self.source))
    self.logger.debug("ROW:\ttbdb_comment: {}".format(self.tbdb_comment))
       
  def complete_row(self):
    """
    This function finishes each row with the rest of the values needed.
    """    
    self.logger.debug("ROW:Within the Row class complete_row function")
    
    if "This mutation is outside the expected region" in self.warning:
      self.logger.debug("ROW:This mutation shouldn't exist! Setting Looker & MDL interpretations of 'NA'")
      self.rationale = "NA"
      self.confidence = "NA"
      self.looker_interpretation = "NA"
      self.mdl_interpretation = "NA"
    
    elif self.who_confidence != "No WHO annotation" and self.who_confidence != "" and self.who_confidence != "NA":
      self.logger.debug("ROW:WHO annotation identified: converting to the appropriate interpretation")
      self.looker_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["looker"]
      self.mdl_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["mdl"]
      self.rationale = "WHO classification"
    
    elif self.who_confidence != "NA":
      self.logger.debug("ROW:No WHO annotation identified: convert with expert rules")
      self.looker_interpretation = self.variant.apply_expert_rules("looker")
      self.mdl_interpretation = self.variant.apply_expert_rules("mdl")
      self.rationale = "Expert rule applied"
      self.confidence = "No WHO annotation"
      
    self.logger.debug("ROW:Interpretation logic applied or skipped; now removing any 'noexpert' suffixes")
    self.describe_rationale()
    self.logger.debug("ROW:Finished completing the row's values, now exiting function")
    
  def rank_annotation(self): 
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
    
  def annotation_to_LIMS(self):
    """
    This function converts the WHO annotation and the target drug
    into returns the LIMS' report file appropriate annotation.
    """
    if self.who_confidence == "Assoc w R":
      return "Mutation(s) associated with resistance to {} detected".format(self.antimicrobial)
    elif (self.who_confidence == "Assoc w R - interim") or (self.who_confidence == "Uncertain significance"):
      return "The detected mutation(s) have uncertain significance. Resistance to {} cannot be ruled out".format(self.antimicrobial)
    # "Not assoc w R" and "Not assoc w R - Interim" and anything else
    else: 
      return "No mutations associated with resistance to {} detected".format(self.antimicrobial)

  def describe_rationale(self):
    """
    This function removes the 'noexpert' suffix in the case where 
    the interpretation logic applied is not considered an expert rule.
    """
    if any(rule in self.looker_interpretation for rule in globals.RULE_TO_RATIONALE.keys()):
      interpretation = self.looker_interpretation[0]
      rule = self.looker_interpretation[1:]
      self.looker_interpretation = interpretation
      self.rationale = globals.RULE_TO_RATIONALE[rule]
      self.confidence = "No WHO annotation"
      
    if any(rule in self.mdl_interpretation for rule in globals.RULE_TO_RATIONALE.keys()):
      interpretation = self.mdl_interpretation[0]
      rule = self.mdl_interpretation[1:]
      self.logger.debug("rule={}".format(rule))
      self.mdl_interpretation = interpretation
      self.logger.debug(globals.RULE_TO_RATIONALE[rule])
      globals.RULE_TO_RATIONALE[rule]
      self.rationale = globals.RULE_TO_RATIONALE[rule]
      self.confidence = "No WHO annotation"
      