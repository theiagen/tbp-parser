from Row import Row
from Variant import Variant
import globals
import pandas as pd
import json

class Laboratorian:
  """
  This class creates the CDPH Laboratorian report.
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
   
  def iterate_section(self, variant_section, row_list):
    for variant in variant_section:
      # create a Variant object and add the origin gene to the GENES_REPORTED set
      variant = Variant(self.logger, variant)
      globals.GENES_REPORTED.add(variant.gene)
      
      # extract all of the annotations for the variant
      variant.extract_annotations()
      
      for annotation_row in variant.annotation_dictionary.values():
        # complete the row objects
        annotation_row.complete_row()
        annotation_row.print()
        
        self.logger.debug("New row: {}".format(annotation_row))
        row_list.append(annotation_row)
        
    return row_list
      
  def create_laboratorian_report(self):
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
    
    self.logger.debug("Creating row_list and the genes_reported list")
    row_list = []
    
    with open(self.input_json) as json_fh:
      input_json = json.load(json_fh)
      
      globals.SAMPLE_NAME = input_json["id"]
      
      row_list = self.iterate_section(input_json["dr_variants"], row_list)
      row_list = self.iterate_section(input_json["other_variants"], row_list)
      
    ### TO-DO: ADD COVERAGE WARNINGS ### 
      
    # add any genes that are missing from the report
    for gene, antimicrobial_drug_names in globals.GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
      for drug_name in antimicrobial_drug_names:
        if gene not in globals.GENES_REPORTED:
          self.logger.debug("Gene {} not in report, adding based on coverage".format(gene))
          row_list.append(Row(self.logger, None, "NA", drug_name, gene))
        else:
          self.logger.debug("Gene {} already in report".format(gene))
          
          ### TO-DO: DOUBLE CHECK THIS LOGIC ###
          
          # check if the gene has at least an R mutation; 
          # applying part 2 of rule 4.2 in the interpretation logic document
          no_r_mutations = set()
          for row in row_list:
            if row.tbprofiler_gene_name == gene:
              if row.looker_interpretation != "R" and "deletion" not in row.warning:
                # if there is no R interpretation for the gene, add it to the no_r_mutations 
                # set, where we will add in the insufficient coverage warning
                no_r_mutations.add(row.tbprofiler_gene_name)
              else:
                # if there is an R interpretation or a deletion for the gene, remove it from the set
                no_r_mutations.discard(row.tbprofiler_gene_name)

          # if the set is not empty, add the insufficient coverage warning
          if len(no_r_mutations) > 0:
            if gene in globals.LOW_DEPTH_OF_COVERAGE_LIST:
              # remove all rows in the row_list entity that belong to this gene
              row_list = [row for row in row_list if row.tbprofiler_gene_name != gene]
              # re-add the row, but with correct information as of 4.2
              row_list.append(Row(self.logger, None, "NA", drug_name, gene))     
     
    self.logger.debug("Creating the dataframe")
    
    # ad row list to DF_LABORATORIAN
    for row in row_list:
      row_dictionary = pd.DataFrame(vars(row), index=[0])
      row_dictionary.drop(["logger", "variant", "who_confidence"], axis=1, inplace=True)
      globals.DF_LABORATORIAN = pd.concat([globals.DF_LABORATORIAN, row_dictionary], ignore_index=True)
      
            
    globals.DF_LABORATORIAN.to_csv("{}.laboratorian_report.csv".format(self.output_prefix), index=False)