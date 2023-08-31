import Variant

"""
This class represents a row in the CDPH Laboratorian report.
"""

class Row(Variant) :
  def __init__(self, variant, samplename, min_depth, coverage_threshold):
    super(variant).__init__()

    # initialize the interpretation fields
    self.looker_interpretation = ""
    self.mdl_interpretation = ""
   
    self.sample_id = samplename
    self.tbprofiler_gene_name = self.gene
    # add gene tier if present in the GENE_TO_TIER look-up dictionary
    if self.gene in globals.GENE_TO_TIER.keys():
      self.gene_tier = globals.GENE_TO_TIER[self.gene]
    else:
      self.gene_tier = "NA"
    
    # fill in the rest of the values needed in the row according to the tbprofiler output contents
    self.tbprofiler_locus_tag = self.locus_tag
    self.tbprofiler_variant_substitution_type = self.type
    self.tbprofiler_variant_substitution_nt = self.nucleotide_change
    self.tbprofiler_variant_substitution_aa = self.protein_change if self.protein_change != "" else "NA"
    self.depth = int(self.depth or 0)
    self.frequency = self.freq
    self.read_support = self.depth*self.frequency
    self.warning = ""
    # Set a warning if the coverage level is below the pre-determined threshold
    # if the mutation is not a deletion, then add the gene to the low coverage list
    if self.depth < int(min_depth) or float(globals.GENE_COVERAGE_DICT[self.gene]) < coverage_threshold:
      self.warning = "Insufficient coverage in locus"
      # I am not sure if this outcome is still wanted
      if "del" in self.nucleotide_change:
        self.warning = "Insufficient coverage in locus (deletion identified)"
      else:
        globals.LOW_DEPTH_OF_COVERAGE_LIST.append(self.gene)


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