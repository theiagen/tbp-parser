
import globals
import pandas as pd
import json
import datetime

class Looker:
  """ 
  This class creates the CDPH Looker report.
  """
  def __init__(self, logger, input_json, output_prefix):
    self.logger = logger
    self.input_json = input_json
    self.output_prefix = output_prefix
  
  def get_lineage_and_id(self):
    """
    Returns the lineage and ID fields for Looker
    """
    with open(self.input_json) as json_fh:
      input_json = json.load(json_fh)
      
      sublineage = input_json["sublin"]
      lineage = "NA"
      if "lineage" in sublineage:
        lineage = sublineage
        ID = "MtBC, not M. bovis"
      elif "BCG" in sublineage:
        ID = "M. bovis BCG"
      elif "bovis" in sublineage:
        ID = "M. bovis, not BCG"
      elif sublineage == "":
        ID = "NA"
      else:
        ID = sublineage
      
      return lineage, ID  
  
  def create_looker_report(self):
    """
    This function recieves the input json and laboratorian report to
    write the Looker report that includes the following information: 
      - sample_id: the sample name
      - for each antimicrobial, indication if resistant (R) or susceptible (S)
    """
    self.logger.info("Within create_looker_report function")
    
    DF_LOOKER = pd.DataFrame({"sample_id": globals.SAMPLE_NAME, "output_seq_method_type": globals.SEQUENCING_METHOD}, index=[0])
    
    # iterate through laboratorian dataframe to extract highest mutation
    for antimicrobial in globals.ANTIMICROBIAL_DRUG_NAME_LIST:
      self.logger.debug("Antimicrobial: {}".format(antimicrobial))
      potential_resistances = globals.DF_LABORATORIAN[globals.DF_LABORATORIAN["antimicrobial"] == antimicrobial]["looker_interpretation"]
      # this is a crazy one liner:
      # basically, it gets the max resistance ranking (R > R-Interim > U > S-Interim > S) for all resistance annotations for a drug
      max_resistance = [annotation for annotation,rank in globals.RESISTANCE_RANKING.items() if rank == max([globals.RESISTANCE_RANKING[interpretation] for interpretation in potential_resistances])]
      DF_LOOKER[antimicrobial] = max_resistance[0]
      self.logger.debug("Resistance: {}".format(max_resistance[0]))
      for gene in globals.ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME[antimicrobial]:
        # indicate warning if any genes failed to achieve 100% coverage_threshold and/or minimum depth  (10x) 
        if DF_LOOKER[antimicrobial][0] != "R" and gene in globals.LOW_DEPTH_OF_COVERAGE_LIST:
          DF_LOOKER[antimicrobial] = "Insufficient coverage for the locus"
            
    # get lineage and ID
    lineage, ID = self.get_lineage_and_id()
    
    # get time stamp
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    
    DF_LOOKER["lineage"] = lineage
    DF_LOOKER["ID"] = ID
    DF_LOOKER["analysis_date"] = current_time
    DF_LOOKER["operator"] = globals.OPERATOR
    
    # write to file
    DF_LOOKER.to_csv("{}.looker_report.csv".format(self.output_prefix), index=False)
        