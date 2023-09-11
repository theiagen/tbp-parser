from Row import Row
from Variant import Variant
import globals
import pandas as pd
import json

class Laboratorian:
  """
  This class creates the CDPH Laboratorian report.
  It has two functions:
    - iterate_section: iterates through each variant in a section and creates
      a Variant in the Variant class for each variant
    - create_laboratorian_report: creates the laboratorian report CSV file
  """
  
  def __init__(self, logger, input_json, output_prefix):
    self.logger = logger
    self.input_json = input_json
    self.output_prefix = output_prefix
   
  def iterate_section(self, variant_section, row_list):
    """
    This function iterate through each variant in a section of the TBProfiler 
    JSON file; for example, it goes through each subsection in the "dr_variants" 
    and "other_variants" sections. It takes each subsection and converts it 
    into an individual Variant object. Then, each annotation within that variant 
    is extracted and converted into a Row object.
    """
    self.logger.info("Within the Laboratorian class iterate_section function")
    
    self.logger.debug("Iterating through the variant section to turn each one into a Variant object")
    for variant in variant_section:
      # create a Variant object and add the origin gene to the global GENES_REPORTED set variable
      variant = Variant(self.logger, variant)
      globals.GENES_REPORTED.add(variant.gene)
      
      # extract all of the annotations for the variant
      variant.extract_annotations()
      
      self.logger.debug("The current variant has {} annotations; now iterating through them".format(len(variant.annotation_dictionary)))
      for annotation_row in variant.annotation_dictionary.values():
        # complete the row objects
        annotation_row.complete_row()
        
        # if in --debug mode, print the annotation row.
        annotation_row.print()
        
        self.logger.debug("New row created! Adding to row_list")
        row_list.append(annotation_row)
        
    self.logger.info("Finished iterating through the variant section, now exiting function")    
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
    self.logger.info("Within the Laboratorian class create_laboratorian_report function")   
    
    row_list = []
    self.logger.debug("Initializing the row_list; contains {} rows".format(len(row_list)))
    
    with open(self.input_json) as json_fh:
      input_json = json.load(json_fh)
      
      globals.SAMPLE_NAME = input_json["id"]
      
      self.logger.debug("About to parse through the variant sections for the sample with name {}".format(globals.SAMPLE_NAME))
      
      row_list = self.iterate_section(input_json["dr_variants"], row_list)
      row_list = self.iterate_section(input_json["other_variants"], row_list)
      
      self.logger.debug("Iteration complete, there are now {} rows".format(len(row_list)))

    self.logger.debug("Now adding any genes that are missing from the report and editing any rows that need to be edited")
    for gene, antimicrobial_drug_names in globals.GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
      for drug_name in antimicrobial_drug_names:
        if gene not in globals.GENES_REPORTED:
          self.logger.debug("Gene {} not in report, now adding it to the report".format(gene))
          row_list.append(Row(self.logger, None, "NA", drug_name, gene))
        else:
          self.logger.debug("Gene {} already in report".format(gene))
          
      if gene in globals.LOW_DEPTH_OF_COVERAGE_LIST:
        self.logger.debug("Checking if the gene ({}) has at least an R mutation if poor coverage, otherwise the row will be overwritten".format(gene))
        no_r_mutations = set()
        r_mutations = set()
        for row in row_list:
          if row.tbprofiler_gene_name == gene:
            # the row contains a mutation, but that mutation is not resistant (whole locus fail point D)
            if (row.looker_interpretation != "R") and any(warning for warning in row.warning if "deletion" in warning) and (row.tbprofiler_variant_substitution_nt != "Insufficient Coverage"):
              no_r_mutations.add(row.tbprofiler_gene_name)
            else:
              r_mutations.add(row.tbprofiler_gene_name)

        # if the no_r_mutations set is not empty but r_mutations is empty, then we want to add an insufficient
        #  coverage warning
        if len(no_r_mutations) > 0 and len(r_mutations) == 0:
            # add a warning to all rows in the row_list entity that belong to this gene
            for row in row_list:
              if row.tbprofiler_gene_name == gene:
                if ("Insufficient coverage in locus (deletion identified)") not in row.warning:
                  row.warning.append("Insufficient coverage in locus")
                  
    self.logger.debug("Creation of rows completed; there are now {} rows".format(len(row_list))) 
        
    # add row list to DF_LABORATORIAN
    for row in row_list:
      row.warning = list(filter(None, row.warning))
      row.warning = ". ".join(row.warning)
      
      # make a temporary dataframe out of the Row object using vars(row) which converts the object into a dictionary
      row_dictionary = pd.DataFrame(vars(row), index=[0])
      row_dictionary.drop(["logger", "variant", "who_confidence"], axis=1, inplace=True)
      globals.DF_LABORATORIAN = pd.concat([globals.DF_LABORATORIAN, row_dictionary], ignore_index=True)
              
    globals.DF_LABORATORIAN.to_csv("{}.laboratorian_report.csv".format(self.output_prefix), index=False)
    self.logger.info("Laboratorian report created, now exiting function")