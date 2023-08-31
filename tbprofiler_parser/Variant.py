import re
import globals
import Annotation

"""
This class represents the variant field from TBProfiler
"""

class Variant:
  def __init__(self, variant=None):
    if variant is not None:
      for key, value in variant.items():
        setattr(self, key, value)
        
      # if the annotation exists, 
      if hasattr(self, "annotation"):
        if len(self.annotation) > 0:
          # turn each annotation into a member of the Annotation class
          annotation_list = []
          for item in self.annotation:
            annotation = Annotation(item["annotation"], item["drug"])
            annotation_list.add(annotation)
          
          # overwrite the annotation field with the list of Annotation objects
          self.annotation = annotation_list

        
      
  
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
      elif "upstream_gene_variant" in self.substitution_type: 
        return "S" if interpretation_destination == "MDL" else "U"
      elif not any(non_ORF in self.nucleotide_change for non_ORF in ["+", "-", "*"]) or self.nucleotide_change.endswith("*"): 
        # if a position includes either +, *, or - it's not in the ORF 
        # UNLESS the * is at the end which means its a premature stop codon
        if self.substitution_type != "synonymous_variant":
          return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
        else:
          return "S"

    elif self.gene == "rrl":
      if (globals.SPECIAL_POSITIONS[self.gene][1][1] <= position_nt <= globals.SPECIAL_POSITIONS[self.gene][1][2]) or (globals.SPECIAL_POSITIONS[self.gene][2][1] <= position_nt <= globals.SPECIAL_POSITIONS[self.gene][2][2]):
        return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
      else:
        return "S" if interpretation_destination == "MDL" else "U"

    # apply expert rules 2.2.1
    elif self.gene in ["katG", "pncA", "ethA", "gid"]: 
      if any(indel_or_stop in self.nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or self.nucleotide_change.endswith("*"):
        return "Uncertain significance" if interpretation_destination == "LIMS" else "U"
      elif (self.substitution_type != "synonymous_variant") or ("upstream_gene_variant" in self.substitution_type):
        return "S" if interpretation_destination == "MDL" else "U"
      else:
        return "S"

    # apply expert rules 2.2.2
    elif self.gene == "rpoB": 
      if globals.SPECIAL_POSITIONS[self.gene][1] <= position_aa <= globals.SPECIAL_POSITIONS[self.gene][2]:
          if self.substitution_type != "synonymous_variant":
            return "Assoc with R" if interpretation_destination == "LIMS" else "R"
          else:
            return "S"   
      elif (self.substitution_type != "synonymous_variant") or ("upstream_gene_variant" in self.substitution_type):
        return "S" if interpretation_destination == "MDL" else "U"
      else:
        return "S"

     # apply rule 3.2 -- not an expert rule
    elif self.gene not in globals.GENE_LIST_OPTION_1 or self.gene not in globals.GENE_LIST_OPTION_2:
      if self.gene == "rrs":
        if position_nt in globals.SPECIAL_POSITIONS[self.gene]:
          return "Unoexpert"
        else:
          return "Snoexpert" if interpretation_destination == "MDL" else "Unoexpert"
      elif (self.substitution_type != "synonymous_variant") or ("upstream_gene_variant" in self.substitution_type):
        return "Snoexpert" if interpretation_destination == "MDL" else "Unoexpert"
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
    