import re
import globals
from Row import Row 

class Variant:
  """
  This class represents a single Variant reported by TBProfiler.
  Note:
  A variant can have either (1) an annotation field. This field could
  (a) have 1 or more annotations, or (b) have none. A variant could also 
  have (2) no annotation field at all.
  
  Within the annotation field, an annotation could have more than one drug
  listed. In addition, it is possible that the same drug can show up twice
  within the same annotation for a single variant.
  
  This class has three functions:
    - extract_annotations: splits up the annotation field into individual
      annotations and creates a Row object for each annotation.
    - apply_expert_rules: applies the expert rules from the CDPH interpretation
      logic document to the variant.
  """
  def __init__(self, logger, variant=None):
    self.logger = logger
    
    # a list containing the various annotations for this variant
    self.annotation_dictionary = {}
    
    if variant is not None:
      for key, value in variant.items():
        setattr(self, key, value)
  
   
  def extract_annotations(self):
    """
    This function takes the annotation field in a variant and splits it into
    its individual parts, creating an Row object for the annotation field. If the
    annotation field is empty or does not exist, the function creates a row
    based off of the gene_associated_drugs field.
    """
    self.logger.debug("Within the Variant class extract_annotations function")
    
    # if possibility 1a (variant has an annotation field with content)
    if hasattr(self, "annotation") and len(self.annotation) > 0:
      self.logger.debug("Starting to turn each annotation into a Row")
      
      # iterate through the annotations
      for item in self.annotation:
        # turn the annotation into a Row object
        annotation = Row(self.logger, self, item["who_confidence"], item["drug"])

        # if this is the first time a drug has been seen, add it to the annotation dictionary
        if annotation.antimicrobial not in self.annotation_dictionary.keys():
          self.logger.debug("This is the first time this drug ({}) has been seen; adding it to the annotation dictionary".format(annotation.antimicrobial))
          self.annotation_dictionary[annotation.antimicrobial] = Row(self.logger, self, annotation.who_confidence, annotation.antimicrobial)
        
        # otherwise, save only the annotation with the more severe WHO confidence (higher value)
        elif annotation.rank_annotation() > self.annotation_dictionary[annotation.antimicrobial].rank_annotation():
          self.logger.debug("A more severe WHO confidence annotation was found for this drug ({}); replacing the old annotation with the new one".format(annotation.antimicrobial))
          self.annotation_dictionary[annotation.antimicrobial] = Row(self.logger, self, annotation.who_confidence, annotation.antimicrobial)
        
      self.logger.debug("The annotation dictionary has been expanded; it now has a length of {}".format(len(self.annotation_dictionary)))
    else:
      # possibilities 1b and 2: the annotation field has no content or the field does not exist
      self.logger.debug("The annotation field has no content or does not exist. Now iterating through gene associated drugs.")
       
      for drug in self.gene_associated_drugs:
        self.annotation_dictionary[drug] = Row(self.logger, self, "No WHO annotation", drug)
  
      self.logger.debug("After iterating through gene_associated_drugs, the annotation dictionary has a length of {}".format(len(self.annotation_dictionary)))

    self.logger.debug("Annotations extracted, now exiting function")
    
  def apply_expert_rules(self, interpretation_destination):
    """
    Apply rules 1-3 from the CDPH interpretation logic document regarding the interpretation of potential resistance mutations.
    """
    self.logger.debug("Within the Variant class apply_expert_rules function")
    
    position_nt = globals.get_position(self.nucleotide_change)
    position_aa = globals.get_position(self.protein_change)

    self.logger.debug("The nucleotide position is {} and the amino acid position is {}".format(position_nt, position_aa))
    
    if self.gene in ["Rv0678", "atpE", "pepQ", "rplC", "mmpL5", "mmpS5"]:        
      self.logger.debug("The gene is {}, now checking if the position requires special consideration".format(self.gene))
     
      # check if position within promoter regions
      if globals.PROMOTER_REGIONS[self.gene][0] <= position_nt <= globals.PROMOTER_REGIONS[self.gene][1]: 
        self.logger.debug("The position is within the promoter region; interpretation is 'U'")
        return "U"
      
      # otherwise, check if it is an upstream gene variant
      elif "upstream_gene_variant" in self.type: 
        self.logger.debug("The position is an upstream gene variant; interpretation is 'S/U'")
        return "S" if interpretation_destination == "mdl" else "U"
      
      # otherwise, check if it is not in the ORF
      elif not any(non_ORF in self.nucleotide_change for non_ORF in ["+", "-", "*"]) or self.nucleotide_change.endswith("*"):
        self.logger.debug("The position is not in the ORF; interpretation is 'S' if it is a synonymous variant or 'U' if it is not")
        # if a position includes either +, *, or - it's not in the ORF, unless 
        #  the * is at the end which means its a premature stop codon
        if self.type != "synonymous_variant":
          return "U"
        else:
          return "S"

    elif self.gene == "rrl":
      self.logger.debug("The gene is rrl, now checking if the position requires special consideration")
      self.logger.debug(globals.SPECIAL_POSITIONS["rrl"][0][1])
      if (globals.SPECIAL_POSITIONS[self.gene][0][0] <= position_nt <= globals.SPECIAL_POSITIONS[self.gene][0][1]) or (globals.SPECIAL_POSITIONS[self.gene][1][0] <= position_nt <= globals.SPECIAL_POSITIONS[self.gene][1][1]):
        self.logger.debug("The position is within the special positions; interpretation is 'U'")
        return "U"
      
      else:
        self.logger.debug("The position is not within the special positions; interpretation is 'S/U'")
        return "S" if interpretation_destination == "mdl" else "U"

    elif self.gene in ["katG", "pncA", "ethA", "gid"]: 
      self.logger.debug("The gene is {}, now checking if the mutation type requires special consideration".format(self.gene))
      
      if any(indel_or_stop in self.nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or self.nucleotide_change.endswith("*"):
        self.logger.debug("The mutation type is an indel, stop, or frameshift codon; interpretation is 'U'")
        return "U"
      
      elif (self.type != "synonymous_variant") or ("upstream_gene_variant" in self.type):
        self.logger.debug("The mutation type is not a synonymous variant or is an upstream gene variant; interpretation is 'S/U'")
        return "S" if interpretation_destination == "mdl" else "U"
      else:
        self.logger.debug("The mutation type is a synonymous variant or is not upstream gene variant; interpretation is 'S'")
        return "S"

    elif self.gene == "rpoB": 
      self.logger.debug("The gene is rpoB, now checking if the position requires special consideration")
      
      if globals.SPECIAL_POSITIONS[self.gene][0] <= position_aa <= globals.SPECIAL_POSITIONS[self.gene][1]:
        self.logger.debug("The position is within the special positions; interpretation is 'R' if nonsynonymous, else 'S'")
        if self.type != "synonymous_variant":
          return "R"
        else:
          return "S"
           
      elif (self.type != "synonymous_variant") or ("upstream_gene_variant" in self.type):
        self.logger.debug("The position is not within the special positions, is not nonsynonymous, or is an upstream gene variant; interpretation is 'S/U'")
        return "S" if interpretation_destination == "mdl" else "U"
      else:
        self.logger.debug("The position is not within the special positions, is synonymous, or is not an upstream gene variant; interpretation is 'S'")
        return "S"

    elif self.gene not in globals.GENE_LIST:
      self.logger.debug("The gene is not in the gene list that requires an expert rule.")
      
      if self.gene == "rrs":
        self.logger.debug("The gene is rrs, now checking if the position requires special consideration")
        
        if position_nt in globals.SPECIAL_POSITIONS[self.gene]:
          self.logger.debug("The position is within the special positions; interpretation is 'U'")
          return "Unoexpert"
        else:
          self.logger.debug("The position is not within the special positions; interpretation is 'S/U'")
          return "Snoexpert" if interpretation_destination == "mdl" else "Unoexpert"
      
      elif (self.type != "synonymous_variant") or ("upstream_gene_variant" in self.type):
        self.logger.debug("The position is not a synonymous variant or is an upstream gene variant; interpretation is 'S/U'")
        return "Snoexpert" if interpretation_destination == "mdl" else "Unoexpert"
      
      else:
        self.logger.debug("The position is a synonymous variant or is not an upstream gene variant; interpretation is 'S'")
        return "Snoexpert"
    return ""