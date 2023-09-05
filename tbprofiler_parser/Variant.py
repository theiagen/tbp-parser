import re
import globals
from Row import Row 

class Variant:
  """
  This class represents the Variant field from TBProfiler.
  Note:
  A variant can have either (1) an annotation field. This field could be
  (b) have 1 or more annotations, or (b) have none. A variant could also 
  have (2) no annotation field at all.
  
  Within the annotation field, an annotation could have more than one drug
  listed. In addition, it is possible that the same drug can show up twice
  within the same annotation for a single variant.
  """
  def __init__(self, logger, variant=None, gene=None, drug_name=None):
    self.logger = logger
    
    # a list containing the various annotations for this variant
    self.annotation_dictionary = {}
    if variant is not None:
      for key, value in variant.items():
        setattr(self, key, value)
   
  def extract_annotations(self):
    """
    This function takes the annotation field in a variant and splits it into
    its individual parts, creating an Annotation class for each part. If the
    annotation field is empty or does not exist, the function creates a row
    based off of the gene_associated_drugs field.
    """
    # if possibility 1a (variant has an annotation field with content)
    if hasattr(self, "annotation") and len(self.annotation) > 0:
      self.logger.debug("Before splitting up the annotations: {}".format(self.annotation))
      
      # turn each annotation into a member of the Annotation class
      for item in self.annotation:
        # if this is the first time a drug has been seen, add it to the annotation dictionary
        annotation = Row(self.logger, self, item["who_confidence"], item["drug"])
        if annotation.antimicrobial not in self.annotation_dictionary.keys():
          self.annotation_dictionary[annotation.antimicrobial] = Row(self.logger, self, annotation.who_confidence, annotation.antimicrobial)
        
        # otherwise, save only the annotaiton with the more severe WHO confidence (higher value)
        elif annotation.rank_annotation() > self.annotation_dictionary[annotation.antimicrobial].rank_annotation():
          self.annotation_dictionary[annotation.antimicrobial] = Row(self.logger, self, annotation.who_confidence, annotation.antimicrobial)
        
      self.logger.debug("After splitting up the annotations: {}".format(self.annotation_dictionary))
    else:
      # possibilities 1b and 2: the annotation field has no content or the field does not exist
      self.logger.debug("No annotations for this variant")
       
      for drug in self.gene_associated_drugs:
        self.annotation_dictionary[drug] = Row(self.logger, self, "No WHO annotation", drug)
  
      self.logger.debug("After iterating through gene_associated_drugs: {}".format(self.annotation_dictionary))
  
  def apply_expert_rules(self, interpretation_destination):
    """
    Apply rules 1-3 from the CDPH Interpretation logic document regarding the interpretation of potential resistance mutations.
    """

    position_nt = self.get_position(self.nucleotide_change)
    position_aa = self.get_position(self.protein_change)

    # apply expert rules 1.2
    if self.gene in ["Rv0678", "atpE", "pepQ", "rplC", "mmpL5", "mmpS5"]:         
      # check if position within promoter regions
      if globals.PROMOTER_REGIONS[self.gene][1] <= position_nt <= globals.PROMOTER_REGIONS[self.gene][2]: 
        return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
      # otherwise, check if it is an upstream gene variant
      elif "upstream_gene_variant" in self.type: 
        return "S" if interpretation_destination == "mdl" else "U"
      elif not any(non_ORF in self.nucleotide_change for non_ORF in ["+", "-", "*"]) or self.nucleotide_change.endswith("*"): 
        # if a position includes either +, *, or - it's not in the ORF 
        # UNLESS the * is at the end which means its a premature stop codon
        if self.type != "synonymous_variant":
          return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
        else:
          return "S"

    elif self.gene == "rrl":
      if (globals.SPECIAL_POSITIONS[self.gene][1][1] <= position_nt <= globals.SPECIAL_POSITIONS[self.gene][1][2]) or (globals.SPECIAL_POSITIONS[self.gene][2][1] <= position_nt <= globals.SPECIAL_POSITIONS[self.gene][2][2]):
        return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
      else:
        return "S" if interpretation_destination == "mdl" else "U"

    # apply expert rules 2.2.1
    elif self.gene in ["katG", "pncA", "ethA", "gid"]: 
      if any(indel_or_stop in self.nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or self.nucleotide_change.endswith("*"):
        return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
      elif (self.type != "synonymous_variant") or ("upstream_gene_variant" in self.type):
        return "S" if interpretation_destination == "mdl" else "U"
      else:
        return "S"

    # apply expert rules 2.2.2
    elif self.gene == "rpoB": 
      if globals.SPECIAL_POSITIONS[self.gene][1] <= position_aa <= globals.SPECIAL_POSITIONS[self.gene][2]:
          if self.type != "synonymous_variant":
            return "Assoc with R" if interpretation_destination == "LIMS" else "R"
          else:
            return "S"   
      elif (self.type != "synonymous_variant") or ("upstream_gene_variant" in self.type):
        return "S" if interpretation_destination == "mdl" else "U"
      else:
        return "S"

     # apply rule 3.2 -- not an expert rule
    elif self.gene not in globals.GENE_LIST_OPTION_1 or self.gene not in globals.GENE_LIST_OPTION_2:
      if self.gene == "rrs":
        if position_nt in globals.SPECIAL_POSITIONS[self.gene]:
          return "Unoexpert"
        else:
          return "Snoexpert" if interpretation_destination == "mdl" else "Unoexpert"
      elif (self.type != "synonymous_variant") or ("upstream_gene_variant" in self.type):
        return "Snoexpert" if interpretation_destination == "mdl" else "Unoexpert"
      else:
        return "Snoexpert"

    return ""

  def get_position(self, mutation):
    """  
    This function recieves a mutation (e.g. 'p.Met291Ile') and
    returns the position (numerical part) as an Int
    """
    pattern = r"\.\D*(\d+)"
    
    match = re.search(pattern, mutation)
    if match:
      return int(match.group(1))
    
    return 0
  
  
  