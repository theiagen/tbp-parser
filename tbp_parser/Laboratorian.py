from Row import Row
from Variant import Variant
import globals as globals_
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
    self.logger.info("LAB:Within the Laboratorian class iterate_section function")
    
    self.logger.debug("LAB:Iterating through the variant section to turn each one into a Variant object")
    for variant in variant_section:
      # create a Variant object and add the origin gene to the global GENES_REPORTED set variable
      variant = Variant(self.logger, variant)
      globals_.GENES_REPORTED.add(variant.gene_name)
      
      # this currently only applies to rpoB -- renaming the gene to the segment name to get the coverage for QC
      if globals_.TNGS and variant.gene_name in globals_.TNGS_REGIONS.keys():
        self.logger.debug("LAB:[tNGS only] checking to see if this is a split primer")
        if len(globals_.TNGS_REGIONS[variant.gene_name]) > 1:
          if isinstance(globals_.TNGS_REGIONS[variant.gene_name], dict):
            
            for segment in globals_.TNGS_REGIONS[variant.gene_name]:
              self.logger.debug("LAB:[tNGS only] checking if variant from {} is found in segment {}".format(variant.gene_name, segment))
              if (globals_.TNGS_REGIONS[variant.gene_name][segment][0] <= variant.pos <= globals_.TNGS_REGIONS[variant.gene_name][segment][1]):
                variant.gene_name_segment = segment
                self.logger.debug("LAB:[tNGS only] variant from {} is found in segment {}; setting gene_name_segment to segment name".format(variant.gene_name, variant.gene_name_segment))
                break
              
              else:
                self.logger.debug("LAB:[tNGS only] variant from {} is NOT found in segment {}".format(variant.gene_name, segment))
            
            if hasattr(variant, "gene_name_segment") is False:
              self.logger.debug("LAB:[tNGS only] This mutation is not in the expected region! This is not within any segments")
              variant.gene_name_segment = "Outside of expected region"
        
      # extract all of the annotations for the variant
      variant.extract_annotations()
      
      self.logger.debug("LAB:The current variant (gene: {}) has {} annotation(s); now iterating through them".format(variant.gene_name, len(variant.annotation_dictionary)))
      for annotation_row in variant.annotation_dictionary.values():
        # complete the row objects
        annotation_row.complete_row()
        
        # if in --debug mode, print the annotation row.
        annotation_row.print()
        
        row_list.append(annotation_row)
        
        # if mutation in mmpS5, mmpL5, or Rv0678, we want to perform further examination but only if the mutation is NOT R based on rule 1.1;
        # otherwise, we want to examine the "aternate_consequences" section and apply rule 1.2 and keep only the highest severity mutation
        # change 2023-12-15: all mmpS/mmpL/mmpR mutations are reported regardless of WHO classification
        if annotation_row.tbprofiler_gene_name in ["mmpS5", "mmpL5", "Rv0678"]:
        # OLD CONDITIONAL included `and !(annotation_row.mdl_interpretation == "R" and annotation_row.rationale == "WHO classification"):`
          row_list = variant.extract_alternate_consequences(annotation_row, row_list)

    self.logger.info("LAB:Finished iterating through the variant section, now exiting function")    
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
      - read_support: the number of reads supporting the mutation (10, depth*frequency)
      - rationale: the rationale for resistance calling (WHO classification, Expert rule)
      - warning: a column reserved for warnings such as low depth of coverage 
      - source: a column used to indicate the resistance source as specified by TBDB
      - tbdb_comment: a column used to include any additional comments as specified by TBDB
    """
    self.logger.info("LAB:Within the Laboratorian class create_laboratorian_report function")   
    
    row_list = []
    self.logger.debug("LAB:Initializing the row_list; contains {} rows".format(len(row_list)))
    
    with open(self.input_json) as json_fh:
      input_json = json.load(json_fh)
      
      globals_.SAMPLE_NAME = input_json["id"]
      
      self.logger.debug("LAB:About to parse through the variant sections for the sample with name {}".format(globals_.SAMPLE_NAME))
      
      row_list = self.iterate_section(input_json["dr_variants"], row_list)
      row_list = self.iterate_section(input_json["other_variants"], row_list)
      
      self.logger.debug("LAB:Iteration complete, there are now {} rows".format(len(row_list)))

    self.logger.debug("LAB:Now adding any genes that are missing from the report and editing any rows that need to be edited")
    for gene, antimicrobial_drug_names in globals_.GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
      for drug_name in antimicrobial_drug_names:
        if gene not in globals_.GENES_REPORTED:
          self.logger.debug("LAB:Gene {} with antimicrobial {} not in report, now adding it to the report".format(gene, drug_name))
          row_list.append(Row(self.logger, None, "NA", drug_name, gene))
        else:
          self.logger.debug("LAB:Gene {} with antimicrobial {} already in report".format(gene, drug_name))
      
      # make a list to add QC fail rows to end of laboratorian report
      reorder_list = [] 

      if gene in globals_.LOW_DEPTH_OF_COVERAGE_LIST:
        self.logger.debug("LAB:Checking if the gene ({}) has poor coverage, if so the row will be overwritten".format(gene))
        
        for row in row_list:       
          if row.tbprofiler_gene_name == gene:
            if row.tbprofiler_gene_name in globals_.GENES_WITH_DELETIONS:
              self.logger.debug("LAB:This gene has a valid deletion so we are removing any locus QC fails as they could be due to the deletion")
              if "Insufficient coverage in locus" in row.warning:
                row.warning.remove("Insufficient coverage in locus")

            if ("Insufficient coverage in locus" in row.warning and 
                "Failed quality in the mutation position" in row.warning and 
                row.mdl_interpretation == "R"):
              self.logger.debug("LAB:This 'R' mutation cannot be trusted due to failing both QC checks; rewriting interpretation to Insufficient Coverage")

              row.mdl_interpretation = "Insufficient Coverage"
              row.looker_interpretation = "Insufficient Coverage"     
              
              reorder_list.append(row)
            
            # whole locus fail point D (non-R mutations with no deletions)
            # add a boolean here to skip R mutations that failed locus coverage if turned on
            elif (row.mdl_interpretation != "R" or globals_.TREAT_R_AS_S):
              if ("del" not in row.tbprofiler_variant_substitution_nt and
                  "Insufficient coverage in locus" in row.warning):
                self.logger.debug("LAB:This mutation is not an 'R' mutation (or is being treated as such) and has bad locus coverage so we are rewriting the interpretation to Insufficient Coverage")

                # overwrite all interpretation values with Insufficient coverage, etc. as per rule 4.2.1.3.2 in the interpretation document
                row.mdl_interpretation = "Insufficient Coverage"
                row.looker_interpretation = "Insufficient Coverage"
                
                if "Insufficient coverage in locus" not in row.warning:
                  row.warning.append("Insufficient coverage in locus")

                reorder_list.append(row)
              else:
                self.logger.debug("LAB:This mutation is not an 'R' mutation (or is being treated as such) but it has decent locus coverage; keeping interpretation as is.")
            else:
              self.logger.debug("LAB:This mutation is an 'R' mutation with decent position quality (or is a deletion) and TREAT_R_AS_S is set to False; keeping interpretation as is.")
          
        # remove rows in reorder_list from row_list and add them to the end of row_list
        for row in reorder_list:
          row_list.remove(row)
          row_list.append(row)
        
    self.logger.debug("LAB:Creation of rows completed; there are now {} rows".format(len(row_list))) 
        
    # add row list to DF_LABORATORIAN
    for row in row_list:
      row.warning = list(filter(None, row.warning))
      row.warning = ". ".join(row.warning)
      
      # make a temporary dataframe out of the Row object using vars(row) which converts the object into a dictionary
      row_dictionary = pd.DataFrame(vars(row), index=[0])
      row_dictionary.drop(["logger", "variant", "who_confidence"], axis=1, inplace=True)
      
      globals_.DF_LABORATORIAN = pd.concat([globals_.DF_LABORATORIAN, row_dictionary], ignore_index=True)
              
    globals_.DF_LABORATORIAN.to_csv("{}.laboratorian_report.csv".format(self.output_prefix), index=False)
    self.logger.info("LAB:Laboratorian report created, now exiting function\n")