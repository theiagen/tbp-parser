import globals

class Row() :
  """
  This class represents a row in the CDPH Laboratorian report.
  """
  def __init__(self, logger, variant, who_confidence, drug, gene_name=None, coverage_threshold=0):
    self.logger = logger
    self.variant = variant
    self.who_confidence = who_confidence
    self.drug = drug
    
    # Initalizing the rest of the columns for the CDPH Laboratorian report
    self.sample_id = ""
    if variant is not None:
      self.tbprofiler_gene_name = self.variant.gene
      self.tbprofiler_locus_tag = self.variant.locus_tag
      self.tbprofiler_variant_substitution_type = self.variant.type
      self.tbprofiler_variant_substitution_nt = self.variant.nucleotide_change
      self.tbprofiler_variant_substitution_aa = self.variant.protein_change
      if self.tbprofiler_variant_substitution_aa == "":
        self.tbprofiler_variant_substitution_aa = "NA"
      self.confidence = self.who_confidence
      self.antimicrobial = self.drug
      self.looker_interpretation = ""
      self.mdl_interpretation = ""
      self.depth = self.variant.depth
      if self.depth is None:
        self.depth = 0
      self.frequency = self.variant.freq
      try:
        self.read_support = self.variant.depth * self.variant.freq
      except:
        self.read_support = 0
      self.rationale = ""
      self.warning = ""
    else:
      self.tbprofiler_gene_name = gene_name
      self.tbprofiler_locus_tag = globals.GENE_TO_LOCUS_TAG[self.tbprofiler_gene_name]
      self.tbprofiler_variant_substitution_nt = "NA"
      self.tbprofiler_variant_substitution_aa = "NA"
      self.confidence = "NA"
      self.antimicrobial = self.drug
      self.depth = "NA"
      self.frequency = "NA"
      self.read_support = "NA"
      self.rationale = "NA"
      self.warning = "NA"
      
      if float(globals.COVERAGE_DICTIONARY[self.tbprofiler_gene_name]) >= coverage_threshold:
        self.tbprofiler_variant_substitution_type = "WT"
        self.tbprofiler_variant_substitution_nt = "WT"
        self.tbprofiler_variant_substitution_aa = "WT"
        self.looker_interpretation = "S"
        self.mdl_interpretation = "WT"
      else:
        self.tbprofiler_variant_substitution_type = "Insufficient Coverage"
        self.looker_interpretation = "Insufficient Coverage"
        self.mdl_interpretation = "Insufficient Coverage"
        self.warning = "Insufficient coverage for the locus"
    
    try:
      self.gene_tier = globals.GENE_TO_TIER[self.tbprofiler_gene_name]
    except:
      self.gene_tier = "NA"
      
  
  def print(self):
    """
    This function prints the row in a readable format.
    """
    self.logger.debug("sample_id: {}".format(self.sample_id))
    self.logger.debug("tbprofiler_gene_name: {}".format(self.tbprofiler_gene_name))
    self.logger.debug("tbprofiler_locus_tag: {}".format(self.tbprofiler_locus_tag))
    self.logger.debug("tbprofiler_variant_substitution_type: {}".format(self.tbprofiler_variant_substitution_type))
    self.logger.debug("tbprofiler_variant_substitution_nt: {}".format(self.tbprofiler_variant_substitution_nt))
    self.logger.debug("tbprofiler_variant_substitution_aa: {}".format(self.tbprofiler_variant_substitution_aa))
    self.logger.debug("confidence: {}".format(self.confidence))
    self.logger.debug("antimicrobial: {}".format(self.antimicrobial))
    self.logger.debug("looker_interpretation: {}".format(self.looker_interpretation))
    self.logger.debug("mdl_interpretation: {}".format(self.mdl_interpretation))
    self.logger.debug("depth: {}".format(self.depth))
    self.logger.debug("frequency: {}".format(self.frequency))
    self.logger.debug("read_support: {}".format(self.read_support))
    self.logger.debug("rationale: {}".format(self.rationale))
    self.logger.debug("warning: {}".format(self.warning))
       
  def complete_row(self, sample_id):
    """
    This function finishes each row with the rest of the values needed.
    """    
    self.logger.info("Within complete_row function")
    
    self.sample_id = sample_id
    self.antimicrobial = self.drug
    
    if self.who_confidence != "No WHO annotation" and self.who_confidence != "" and self.who_confidence != "NA":
      self.logger.info("WHO annotation identified: convert to interpretation logic")
      self.looker_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["looker"]
      self.mdl_interpretation = globals.ANNOTATION_TO_INTERPRETATION[self.who_confidence]["mdl"]
      self.rationale = "WHO classification"
    
    elif self.who_confidence != "NA":
      self.logger.info("No WHO annotation identified: convert with expert rules")
      self.looker_interpretation = self.variant.apply_expert_rules("looker")
      self.mdl_interpretation = self.variant.apply_expert_rules("mdl")
      self.rationale = "Expert rule applied"
      
    self.logger.info("Interpretation logic applied or skipped")
    
    self.remove_no_expert()
    
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
      return "Mutation(s) associated with resistance to {} detected".format(self.drug)
    elif (self.who_confidence == "Assoc w R - interim") or (self.who_confidence == "Uncertain significance"):
      return "The detected mutation(s) have uncertain significance. Resistance to {} cannot be ruled out".format(self.drug)
    # "Not assoc w R" and "Not assoc w R - Interim" and anything else
    else: 
      return "No mutations associated with resistance to {} detected".format(self.drug)

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