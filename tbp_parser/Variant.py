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
  
    if hasattr(self, "gene") and self.gene_name == "mmpR5":
      self.logger.debug("VAR:The gene is mmpR5, renaming to Rv0678")
      self.gene_name = "Rv0678"
  
  def extract_alternate_consequences(self, parent_row, row_list):
    """
    This function takes apart the alternate_consequences field for variants
    in the mmpS5, mmpL5, and Rv0678 genes. It creates a Row object for each, 
    and then determines which one has the highest severity. The highest
    severity mutation is returned.
    """
    self.logger.debug("VAR:Within the Variant class extract_alternate_consequences function")

    # check if alternate_consequences section exists, and if so, enter it
    if hasattr(self, "consequences") and len(self.consequences) > 0:
      self.logger.debug("VAR:Now iterating through alternate_consequences")
      
      for item in self.consequences:
        var = Variant(self.logger, item)
        if var.gene_name == "mmpR5":
          var.gene_name = "Rv0678"
        
        if var.gene_name != parent_row.tbprofiler_gene_name:
          globals.GENES_REPORTED.add(var.gene_name)
          self.logger.debug("VAR:This alternate consequence is not from the same gene; creating a new row")
        
          row = Row(self.logger, var, parent_row.who_confidence, parent_row.antimicrobial, var.gene_name, parent_row.depth, parent_row.frequency)
          row.source = parent_row.source
          row.tbdb_comment = parent_row.tbdb_comment
          row.complete_row()
      
          row_list.append(row)
        else:
          self.logger.debug("VAR:This alternate consequence is from the same gene; we do not want to create a new row")
          continue
           
    return row_list
   
  def extract_annotations(self):
    """
    This function takes the annotation field in a variant and splits it into
    its individual parts, creating an Row object for the annotation field. If the
    annotation field is empty or does not exist, the function creates a row
    based off of the gene_associated_drugs field.
    """
    self.logger.debug("VAR:Within the Variant class extract_annotations function")
    
    # if possibility 1a (variant has an annotation field with content)
    if hasattr(self, "annotation") and len(self.annotation) > 0:
      self.logger.debug("VAR:Starting to turn each annotation into a Row")
      
      # create a list of the drugs associated with the gene to check if all drugs are reported
      gene_associated_drug_list = self.gene_associated_drugs
      gene_associated_drug_list = ["rifampin" if drug == "rifampicin" else drug for drug in gene_associated_drug_list] # rename rifampicin to rifampin

      # iterate through the annotations
      for item in self.annotation:      
        # turn the annotation into a Row object
        annotation = Row(self.logger, self, item["confidence"], item["drug"], source=item["source"], tbdb_comment=item["comment"])
        
        # if this is the first time a drug has been seen, add it to the annotation dictionary
        if annotation.antimicrobial not in self.annotation_dictionary.keys():
          self.logger.debug("VAR:This is the first time this drug ({}) has been seen; adding it to the annotation dictionary".format(annotation.antimicrobial))
          self.annotation_dictionary[annotation.antimicrobial] = Row(self.logger, self, annotation.confidence, annotation.antimicrobial, source=annotation.source, tbdb_comment=annotation.tbdb_comment)

          # if the drug is in the gene associated drug list, remove it because it has been accounted for
          if annotation.antimicrobial in gene_associated_drug_list:
            gene_associated_drug_list.remove(annotation.antimicrobial)
            self.logger.debug("VAR:The gene associated drug list: {}".format(gene_associated_drug_list))
        
        # otherwise, save only the annotation with the more severe WHO confidence (higher value)
        elif annotation.rank_annotation() > self.annotation_dictionary[annotation.antimicrobial].rank_annotation():
          self.logger.debug("VAR:A more severe WHO confidence annotation was found for this drug ({}); replacing the old annotation with the new one".format(annotation.antimicrobial))
          self.annotation_dictionary[annotation.antimicrobial] = Row(self.logger, self, annotation.confidence, annotation.antimicrobial, source=annotation.source, tbdb_comment=annotation.tbdb_comment)
        
      self.logger.debug("VAR:The annotation dictionary has been expanded; it now has a length of {}".format(len(self.annotation_dictionary)))

      # add any missing drugs from the annotation list that are found on the gene associated drug list to the annotation dictionary 
      for drug in gene_associated_drug_list:  
        self.logger.debug("VAR:The drug ({}) was not found in the annotation dictionary; adding it with a WHO confidence of 'No WHO annotation'".format(drug))
        self.annotation_dictionary[drug] = Row(self.logger, self, "No WHO annotation", drug)
        
      self.logger.debug("VAR:The annotation dictionary has all gene associated drugs included; it now has a length of {}".format(len(self.annotation_dictionary)))
    else:
      # possibilities 1b and 2: the annotation field has no content or the field does not exist
      self.logger.debug("VAR:The annotation field has no content or does not exist. Now iterating through gene associated drugs.")
       
      for drug in self.gene_associated_drugs:
        self.annotation_dictionary[drug] = Row(self.logger, self, "No WHO annotation", drug)
  
      self.logger.debug("VAR:After iterating through gene_associated_drugs, the annotation dictionary has a length of {}".format(len(self.annotation_dictionary)))

    self.logger.debug("VAR:Annotations extracted, now exiting function")
    
  def apply_expert_rules(self, interpretation_destination):
    """
    Apply rules 1-3 from the CDPH interpretation logic document regarding the interpretation of potential resistance mutations.
    """
    self.logger.debug("VAR:Within the Variant class apply_expert_rules function")
    
    position_nt = globals.get_position(self.nucleotide_change)
    position_aa = globals.get_position(self.protein_change)

    self.logger.debug("VAR:The nucleotide position is {} and the amino acid position is {}".format(position_nt, position_aa))
    
    if self.gene_name in ["Rv0678", "atpE", "pepQ", "rplC", "mmpL5", "mmpS5"]:        
      self.logger.debug("VAR:The gene is {}, now checking if the position requires special consideration under rule 1.2".format(self.gene_name))
     
      # check if position within promoter regions
      if self.gene_name in globals.PROMOTER_REGIONS.keys():
        if ((len(position_nt) > 1 and 
             (any([x in range(globals.PROMOTER_REGIONS[self.gene_name][0], globals.PROMOTER_REGIONS[self.gene_name][1]) for x in position_nt]) 
              or any([x in range(position_nt[0], position_nt[1]) for x in globals.PROMOTER_REGIONS[self.gene_name]]))) 
            or (globals.PROMOTER_REGIONS[self.gene_name][0] <= position_nt[0] <= globals.PROMOTER_REGIONS[self.gene_name][1])):
          self.logger.debug("VAR:The position is within the promoter region; interpretation is 'U'")
          return "Urule1.2"
      
      # otherwise, check if it is an upstream gene variant
      if "upstream_gene_variant" in self.type: 
        self.logger.debug("VAR:The position is an upstream gene variant; interpretation is 'S/U'")
        return "Srule1.2" if interpretation_destination == "mdl" else "Urule1.2"
      
      # otherwise, check if it is not in the ORF
      if ((not any(non_ORF in self.nucleotide_change for non_ORF in ["+", "-", "*"]) or self.nucleotide_change.endswith("*")) 
          or (not any(non_ORF in self.protein_change for non_ORF in ["+", "-", "*"]) or self.protein_change.endswith("*"))):
        self.logger.debug("VAR:The position is not in the ORF; interpretation is 'S' if it is a synonymous variant or 'U' if it is not")
        # if a position includes either +, *, or - it's not in the ORF, unless 
        #  the * is at the end which means its a premature stop codon
        if self.type != "synonymous_variant":
          return "Urule1.2"
        else:
          return "Srule1.2"

    elif self.gene_name == "rrl":
      self.logger.debug("VAR:The gene is rrl, now checking if the position requires special consideration")
      
      # there is surely a better way to do this, but this works. essentially, I'm checking if the position(s) of the rrl mutation are within the special position range
      if ((len(position_nt) > 1 and 
           (any([x in range(globals.SPECIAL_POSITIONS[self.gene_name][0][0], globals.SPECIAL_POSITIONS[self.gene_name][0][1]) for x in position_nt]) 
            or any([x in range(globals.SPECIAL_POSITIONS[self.gene_name][1][0], globals.SPECIAL_POSITIONS[self.gene_name][1][1]) for x in position_nt]) 
            or any([x in range(position_nt[0], position_nt[1]) for x in globals.SPECIAL_POSITIONS[self.gene_name][0]]) 
            or any([x in range(position_nt[0], position_nt[1]) for x in globals.SPECIAL_POSITIONS[self.gene_name][1]]))) 
          or ((globals.SPECIAL_POSITIONS[self.gene_name][0][0] <= position_nt[0] <= globals.SPECIAL_POSITIONS[self.gene_name][0][1]) 
              or (globals.SPECIAL_POSITIONS[self.gene_name][1][0] <= position_nt[0] <= globals.SPECIAL_POSITIONS[self.gene_name][1][1]))):
        self.logger.debug("VAR:The position is within the special positions; interpretation is 'U'")
        return "Urule1.2"
      
      else:
        self.logger.debug("VAR:The position is not within the special positions; interpretation is 'S/U'")
        return "Srule1.2" if interpretation_destination == "mdl" else "Urule1.2"

    elif self.gene_name in ["katG", "pncA", "ethA", "gid"]: 
      self.logger.debug("VAR:The gene is {}, now checking if the mutation type requires special consideration under rule 2.2".format(self.gene_name))

      if ((any(indel_or_stop in self.nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) 
           or self.nucleotide_change.endswith("*")) 
          or (any(indel_or_stop in self.protein_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) 
              or self.protein_change.endswith("*"))):
        if any([int(position) > -30 for position in position_nt]): 
          self.logger.debug("VAR:The mutation type is an indel, stop, or frameshift codon and within 30 nt of the start codon; interpretation is 'R'")
          return "Rrule2.2.1"
        else:
          self.logger.debug("VAR:The mutation type is an indel, stop, or frameshift codon but is past 30 nt of the start codon; interpretation depends on variant type")
      
      if (self.type == "synonymous_variant"):
        self.logger.debug("VAR:The mutation type is a synonymous variant; interpretation is 'S'")
        return "Srule2.2.1"
      elif ("upstream_gene_variant" in self.type):
        self.logger.debug("VAR:The mutation type is a NONsynonymous variant, and is an upstream gene variant; interpretation is 'S/U'")
        return "Srule2.2.1" if interpretation_destination == "mdl" else "Urule2.2.1"
      else:
        self.logger.debug("VAR:The mutation type IS a NONsynonymous variant, and is NOT an upstream gene variant; interpretation is 'U'")
        return "Urule2.2.1"

    # rules 2.2.2.1 and 3.2.2 & 3.2.3
    elif self.gene_name in ["gyrA", "gyrB", "rpoB"]: 
      self.logger.debug("VAR:The gene is {}, now checking if the position requires special consideration".format(self.gene_name))
      self.logger.debug("VAR:length of aa: {}".format(len(position_aa)))
      self.logger.debug("VAR:SEPCIAL POSITIONS: {}".format(globals.SPECIAL_POSITIONS[self.gene_name]))
      
      # explanation of following giant conditional:
      # RRDR belongs to rpoB, QRDR belongs to gyrA and gyrB
      # if there are more than 2 AA positions, it's likely a deletion or frameshift thing. 
      #  if that is the case, we want to check to see if any of those positions are within (R/Q)RDR *or* if (R/Q)RDR falls within those positions
      #  to see if within (R/Q)RDR, we check if any of the AA positions are within the (R/Q)RDR range
      #  to see if (R/Q)RDR falls within, we check (R/Q)RDR start/end positions are within the AA position range
      # otherwise, check if the single AA position is within (R/Q)RDR positions
      if ((len(position_aa) > 1 and 
           (any([x in range(globals.SPECIAL_POSITIONS[self.gene_name][0], globals.SPECIAL_POSITIONS[self.gene_name][1]) for x in position_aa]) 
            or any([x in range(position_aa[0], position_aa[1]) for x in globals.SPECIAL_POSITIONS[self.gene_name]]))) 
          or ((globals.SPECIAL_POSITIONS[self.gene_name][0] <= position_aa[0] <= globals.SPECIAL_POSITIONS[self.gene_name][1]))):
        self.logger.debug("VAR:The position is within the special positions; interpretation is 'R' if rpoB (or 'U' if not) and nonsynonymous, else 'S'")
        
        if self.gene_name == "rpoB":
          if self.type == "synonymous_variant":
            return "Srule2.2.2.1" 
          else:
            return "Rrule2.2.2.1"
          
        if self.gene_name == "gyrA":
          if self.type != "synonymous_variant":
            return "Urule3.2.2"
          
        if self.gene_name == "gyrB":
          if self.type != "synonymous_variant":
            return "Urule3.2.3"
      
      elif (self.type == "synonymous_variant"):
        self.logger.debug("VAR:The position is not within the special positions and is synonymous; interpretation is 'S'")
        return "Srule2.2.2.2" if self.gene_name == "rpoB" else "Srule3.2.4"
      elif ("upstream_gene_variant" in self.type):
        self.logger.debug("VAR:The position is not within the special positions, is nonsynomyous but is an upstream gene variant; interpretation is 'S/U'")
        if self.gene_name == "rpoB":
          return "Srule2.2.2.2" if interpretation_destination == "mdl" else "Urule2.2.2.2"
        else:
          return "Srule3.2.4" if interpretation_destination == "mdl" else "Urule3.2.4"
      else:
        self.logger.debug("VAR:The position is not within the special positions, is nonsynonymous and is NOT an upstream gene variant; interpretation is 'U'")
        return "Urule2.2.2.2" if self.gene_name == "rpoB" else "Urule3.2.4"

    elif self.gene_name not in globals.GENE_LIST:
      self.logger.debug("VAR:The gene is not in the gene list that requires an expert rule.")
      
      if self.gene_name == "rrs":
        self.logger.debug("VAR:The gene is rrs, now checking if the position requires special consideration")
        
        if any(map(lambda position: position in position_nt, globals.SPECIAL_POSITIONS[self.gene_name])):
          self.logger.debug("VAR:The position is within the special positions; interpretation is 'U'")
          return "Urule3.2.1"
        else:
          self.logger.debug("VAR:The position is not within the special positions; interpretation is 'S/U'")
          return "Srule3.2.1" if interpretation_destination == "mdl" else "Urule3.2.1"
      
      # rule 3.2.4: all remaining scenarios not covered above
      elif (self.type == "synonymous_variant"):
        self.logger.debug("VAR:The mutation is synonymous; interpretation is 'S'")
        return "Srule3.2.4"
      elif ("upstream_gene_variant"  in self.type):
        self.logger.debug("VAR:The mutation is a nonsynonymous variant but is an upstream gene variant; interpretation is 'S/U'")
        return "Srule3.2.4" if interpretation_destination == "mdl" else "Urule3.2.4"
      else:
        self.logger.debug("VAR:The mutation is a nonsynonymous variant or is an upstream gene variant; interpretation is 'U'")
        return "Urule3.2.4"
      
    return ""