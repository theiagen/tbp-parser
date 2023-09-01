from Row import Row
from Variant import Variant
import pandas as pd
import json

class Laboratorian:
  """
  This class creates the laboratorian report.
  """
  def __init__(self, logger, input_json, output_prefix):
    self.logger = logger
    self.input_json = input_json
    self.output_prefix = output_prefix

  def separate_variants(self, variant_section):
    """
    This function separates each variant in a section into its own
    individual Variant class
    """
    self.logger.info("Within separate_variants function")
    for variant in variant_section:
      variant = Variant(variant)
      #globals.VARIANTS.append(variant)
    

  def create_laboratorian_report(self, min_depth, coverage_threshold):
    """
    This function creates the laboratorian report, which is a CSV file
    containing the following information for each mutation in the input JSON:
      - sample_id: the sample name
      - tbprofiler_gene_name: the gene name
      - tbprofiler_locus_tag: the locus tag
      - tbprofiler_variant_substitution_type: the variant substitution type (missense_variant, upstream_gene_variant...)
      - tbprofiler_variant_substitution_nt: the nucleotide substitution (c.1349C>G)
      - tbprofiler_variant_substitution_aa: the amino acid substitution (p.Ser450Trp)
      - confidence: the tbprofiler annotation regarding resistance (Not assoc w R, Uncertain significance...)
      - antimicrobial: the antimicrobial drug the mutation confers resistance to (streptomycin, rifampin...)
      - looker_interpretation: the interpretation of resistance for the CDPH Looker report (R, R-interim, U, S, S-interim)
      - mdl_interpretation: the MDL interpretation of resistance (R, S, U)
      - depth: the depth of coverage at the mutation site (100)
      - frequency: the frequency of mutation at the site (1)
      - read_support: the number of reads supporting the mutation (100, depth*frequency)
      - rationale: the rationale for resistance calling (WHO classification, Expert rule)
      - warning: a column reserved for warnings such as low depth of coverage 
    """
    self.logger.info("Within create_laboratorian_report function")
    
    df_laboratorian = pd.DataFrame(columns = [
      "sample_id", "tbprofiler_gene_name", "tbprofiler_locus_tag", "tbprofiler_variant_substitution_type", 
      "tbprofiler_variant_substitution_nt", "tbprofiler_variant_substitution_aa", "confidence", "antimicrobial",
      "looker_interpretation", "mdl_interpretation", "depth", "frequency", "read_support", "rationale", "warning"
    ])
    
    self.logger.debug("Creating row_list and the genes_reported list")
    row_list = []
    
    with open(self.input_json) as json_fh:
      input_json = json.load(json_fh)
      
      samplename = input_json["id"]
      
      for variant in input_json["dr_variants"]:
        # create a Variant object
        variant = Variant(self.logger, variant)
        variant.extract_annotations()
        
        for annotation_row in variant.annotation_dictionary.values():
          # complete the row objectsSW
          annotation_row.complete_row()
          
          self.logger.debug("New row: {}".format(annotation_row))
          row_list.append(annotation_row)
        # create a row object
        #row = Row(self.logger, variant, samplename, min_depth)
        #row.complete_row()
        
        #self.logger.debug("New row: {}".format(row))
        #row_list.append(row)
    
    self.logger.debug("Creating the dataframe")
    # convert row list into a dataframe
    
    df_laboratorian = df_laboratorian.add(row_list, ignore_index=True)
    df_laboratorian.to_csv("{}.laboratorian_report.csv".format(self.output_prefix), index=False)