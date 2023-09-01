import tbprofiler_parser.Variant as Variant
import tbprofiler_parser.globals as globals
class Row() :
  """
  This class represents a row in the CDPH Laboratorian report.
  """
  def __init__(self, logger, variant, who_confidence, drug):
    self.logger = logger
    self.variant = variant
    self.who_confidence = who_confidence
    self.drug = drug
    
    # Initalizing the rest of the columns for the CDPH Laboratorian report
    self.sample_id = ""
    self.tbprofiler_gene_name = self.variant.gene
    self.tbprofiler_locus_tag = self.variant.locus_tag
    self.tbprofiler_variant_substitution_type = self.variant.substitution_type
    self.tbprofiler_variant_substitution_nt = self.variant.substitution_nt
    self.tbprofiler_variant_substitution_aa = self.variant.substitution_aa
    self.confidence = ""
    self.antimicrobial = self.drug
    self.looker_interpretation = ""
    self.mdl_interpretation = ""
    self.depth = self.variant.depth
    self.frequency = self.variant.frequency
    self.read_support = self.variant.read_support
    self.rationale = ""
    self.warning = ""
   
    
  def complete_row(self):
    """
    This function finishes each row with the rest of the values needed.
    """
    self.logger.info("Within complete_row function")
    
    self.confidence = "No WHO annotation" if self.variant.who_confidence == "No WHO annotation" else self.variant.who_confidence
    self.antimicrobial = self.drug
    
    if self.who_confidence != "No WHO annotation":
      self.logger.info("WHO annotation identified: convert to interpretation logic")
      self.looker_interpretation = globals.ANNOTATION_TO_INTERPRETATION(self.who_confidence, "looker")
      self.mdl_interpretation = globals.ANNOTATION_TO_INTERPRETATION(self.who_confidence, "mdl")
      self.rationale = "WHO classificaiton"
    
    else:
      self.logger.info("No WHO annotation identified: convert with expert rules")
      self.looker_interpretation = self.variant.apply_expert_rules("looker")
      self.mdl_interpretation = self.variant.apply_expert_rules("mdl")
      self.rationale = "Expert rule applied"
      
    self.logger.info("Interpretation logic applied")
    
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
    
    if "noexpert" in self.mdl_interpretation:
      interpretation = self.mdl_interpretation
      self.mdl_interpretation = interpretation.replace("noexpert", "")
      self.rationale = "No WHO annotation or expert rule"