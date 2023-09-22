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
    - remove_no_expert: removes the 'noexpert' suffix in the case where
      the interpretation logic applied is not considered an expert rule
  """
  
  def __init__(self, logger, variant, who_confidence, drug, gene_name=None):
    self.logger = logger
    self.logger.debug("Within the Row class __init__ function")
    
    self.variant = variant
    self.who_confidence = who_confidence
    self.antimicrobial = drug
  
    if gene_name != "test":
      self.antimicrobial = self.antimicrobial.replace("rifampicin", "rifampin")

      self.sample_id = globals.SAMPLE_NAME

      # Initalizing the rest of the columns for the CDPH Laboratorian report
      # for when the variant is in the JSON file
      if variant is not None:
        self.logger.debug("Initalizing the Row object, the variant has been supplied.")
        self.tbprofiler_gene_name = self.variant.gene
        self.tbprofiler_locus_tag = self.variant.locus_tag
        self.tbprofiler_variant_substitution_type = self.variant.type
        self.tbprofiler_variant_substitution_nt = self.variant.nucleotide_change
        self.tbprofiler_variant_substitution_aa = self.variant.protein_change
        # change blank aa substitutions to NA
        if self.tbprofiler_variant_substitution_aa == "":
          self.tbprofiler_variant_substitution_aa = "NA"
        self.confidence = self.who_confidence
        self.looker_interpretation = ""
        self.mdl_interpretation = ""
        self.depth = int(self.variant.depth or 0)
        self.frequency = self.variant.freq
        # avoid division by zero errors
        try:
          self.read_support = self.variant.depth * self.variant.freq
        except:
          self.read_support = 0
        self.rationale = ""
        self.warning = []
        
        if (self.depth < globals.MIN_DEPTH) or (float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) < globals.COVERAGE_THRESHOLD):
          self.logger.debug("The depth of coverage for this variant is {} and the coverage for the gene is {}; applying a warning".format(self.depth, globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]))
          if (float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) < globals.COVERAGE_THRESHOLD):
            globals.LOW_DEPTH_OF_COVERAGE_LIST.append(self.tbprofiler_gene_name)
            
            if "del" in self.tbprofiler_variant_substitution_nt or self.tbprofiler_gene_name in globals.GENES_WITH_DELETIONS:
              self.warning.append("Insufficient coverage in locus (deletion identified)")
              globals.GENES_WITH_DELETIONS.add(self.tbprofiler_gene_name)
            else:
              self.warning.append("Insufficient coverage in locus")
         
          if (self.depth < globals.MIN_DEPTH or float(self.frequency) < 0.10 or self.read_support < 10) and "del" not in self.tbprofiler_variant_substitution_nt:
            globals.MUTATION_FAIL_LIST.append(self.tbprofiler_variant_substitution_nt)
            self.warning.append("Failed quality in the mutation position")
        else:
          self.warning = [""]
      # otherwise, the variant does not appear in the JSON file and default NA/WT values
      # need to be supplied
      else:
        self.logger.debug("Initializing the Row object, the variant has no information supplied. Defaulting to NA or WT values.")
        if gene_name == "mmpR5":
          self.tbprofiler_gene_name = "Rv0678"
        else: 
          self.tbprofiler_gene_name = gene_name
        self.tbprofiler_locus_tag = globals.GENE_TO_LOCUS_TAG[self.tbprofiler_gene_name]
        self.tbprofiler_variant_substitution_nt = "NA"
        self.tbprofiler_variant_substitution_aa = "NA"
        self.confidence = "NA"
        self.depth = "NA"
        self.frequency = "NA"
        self.read_support = "NA"
        self.rationale = "NA"
        self.warning = [""]
        
        # check to see if we need to apply a coverage warning (whole locus fail point c)
        if float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) >= globals.COVERAGE_THRESHOLD:
          self.logger.debug("A warning does not need to be applied; setting the variant information to WT")
          self.tbprofiler_variant_substitution_type = "WT"
          self.tbprofiler_variant_substitution_nt = "WT"
          self.tbprofiler_variant_substitution_aa = "WT"
          self.looker_interpretation = "S"
          self.mdl_interpretation = "WT"
        else:
          self.logger.debug("A warning needs to be applied; setting the variant information to insufficient coverage")
          self.tbprofiler_variant_substitution_type = "Insufficient Coverage"
          self.looker_interpretation = "Insufficient Coverage"
          self.mdl_interpretation = "Insufficient Coverage"
          self.warning.append("Insufficient coverage in locus")
      
      try:
        self.gene_tier = globals.GENE_TO_TIER[self.tbprofiler_gene_name]
      except:
        self.gene_tier = "NA"
      
      self.logger.debug("Row object initialized, exiting __init__ function")
      
  def print(self):
    """
    This function prints the row in a readable format.
    """
    self.logger.debug("Now printing the row in a readable format:")
    self.logger.debug("\tsample_id: {}".format(self.sample_id))
    self.logger.debug("\ttbprofiler_gene_name: {}".format(self.tbprofiler_gene_name))
    self.logger.debug("\ttbprofiler_locus_tag: {}".format(self.tbprofiler_locus_tag))
    self.logger.debug("\ttbprofiler_variant_substitution_type: {}".format(self.tbprofiler_variant_substitution_type))
    self.logger.debug("\ttbprofiler_variant_substitution_nt: {}".format(self.tbprofiler_variant_substitution_nt))
    self.logger.debug("\ttbprofiler_variant_substitution_aa: {}".format(self.tbprofiler_variant_substitution_aa))
    self.logger.debug("\tconfidence: {}".format(self.confidence))
    self.logger.debug("\tantimicrobial: {}".format(self.antimicrobial))
    self.logger.debug("\tlooker_interpretation: {}".format(self.looker_interpretation))
    self.logger.debug("\tmdl_interpretation: {}".format(self.mdl_interpretation))
    self.logger.debug("\tdepth: {}".format(self.depth))
    self.logger.debug("\tfrequency: {}".format(self.frequency))
    self.logger.debug("\tread_support: {}".format(self.read_support))
    self.logger.debug("\trationale: {}".format(self.rationale))
    self.logger.debug("\twarning: {}".format(self.warning))
       
  def complete_row(self):
    """
    This function finishes each row with the rest of the values needed.
    """    
    self.logger.debug("Within the Row class complete_row function")
    
    if self.who_confidence != "No WHO annotation" and self.who_confidence != "" and self.who_confidence != "NA":
      self.logger.debug("WHO annotation identified: converting to the appropriate interpretation")
      self.looker_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["looker"]
      if self.tbprofiler_gene_name in globals.GENE_LIST_MDL:
        self.mdl_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["mdl-ingenelist1"]
      else:
        self.mdl_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["mdl"]
      self.rationale = "WHO classification"
    
    elif self.who_confidence != "NA":
      self.logger.debug("No WHO annotation identified: convert with expert rules")
      self.looker_interpretation = self.variant.apply_expert_rules("looker")
      self.mdl_interpretation = self.variant.apply_expert_rules("mdl")
      self.rationale = "Expert rule applied"
      
    self.logger.debug("Interpretation logic applied or skipped; now removing any 'noexpert' suffixes")
    self.remove_no_expert()
    
    self.logger.debug("Finished completing the row's values, now exiting function")
    
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

  def remove_no_expert(self):
    """
    This function removes the 'noexpert' suffix in the case where 
    the interpretation logic applied is not considered an expert rule.
    """
    if "noexpert" in self.looker_interpretation:
      interpretation = self.looker_interpretation
      self.looker_interpretation = interpretation.replace("noexpert", "")
      self.rationale = "No WHO annotation or expert rule"
      self.confidence = "No WHO annotation"
      
    if "noexpert" in self.mdl_interpretation:
      interpretation = self.mdl_interpretation
      self.mdl_interpretation = interpretation.replace("noexpert", "")
      self.rationale = "No WHO annotation or expert rule"
      self.confidence = "No WHO annotation"