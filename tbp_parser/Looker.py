
import globals as globals_
import pandas as pd
import datetime

class Looker:
    """ 
    This class creates the CDPH Looker report.

    It has one function:
        - create_looker_report: creates the Looker report CSV file 
    """
    
    RESISTANCE_RANKING = {
        "R": 6,
        "R-Interim": 5,
        "U": 4,
        "S-Interim": 3,
        "S": 2,
        "WT": 1,
        "Insufficient Coverage": 0,
        "NA": -1 # outside the expected region
    }    
    """A dictionary ranking the resistance annotations from highest to lowest priority."""
    
    
    def __init__(self, logger, output_prefix, DF_LABORATORIAN, LOW_DEPTH_OF_COVERAGE_LIST, GENES_WITH_VALID_DELETIONS, GENE_TO_ANTIMICROBIAL_DRUG_NAME):
        self.logger = logger
        self.output_prefix = output_prefix
        self.DF_LABORATORIAN = DF_LABORATORIAN
        self.LOW_DEPTH_OF_COVERAGE_LIST = LOW_DEPTH_OF_COVERAGE_LIST
        self.GENES_WITH_VALID_DELETIONS = GENES_WITH_VALID_DELETIONS
        self.GENE_TO_ANTIMICROBIAL_DRUG_NAME = GENE_TO_ANTIMICROBIAL_DRUG_NAME

    def create_looker_report(self, SAMPLE_NAME, SEQUENCING_METHOD, LINEAGE, LINEAGE_ENGLISH, OPERATOR) -> None:
        """
        This function recieves the input json and laboratorian report to
        write the Looker report that includes the following information: 
            - sample_id: the sample name
            - for each antimicrobial, indication if resistant (R) or susceptible (S)
        """
        DF_LOOKER = pd.DataFrame({
            "sample_id": SAMPLE_NAME, 
          "output_seq_method_type": SEQUENCING_METHOD
            }, index=[0])

        # reverse the GENE_TO_ANTIMICROBIAL_DRUG_NAME dictionary to get antimicrobial drug names
        drugs_to_genes = {}
        
        for gene, drugs in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME.items():
            for drug in drugs:
                if drug not in drugs_to_genes:
                    drugs_to_genes[drug] = []
                drugs_to_genes[drug].append(gene)
                

        # iterate through laboratorian dataframe to extract highest mutation
        for antimicrobial in drugs_to_genes.keys():
            potential_looker_resistances = self.DF_LABORATORIAN[self.DF_LABORATORIAN["antimicrobial"] == antimicrobial]["looker_interpretation"]

            try:
                max_looker_resistance = max(potential_looker_resistances, key=lambda x: self.RESISTANCE_RANKING[x])
            except:
                max_looker_resistance = "NA"
                
            DF_LOOKER[antimicrobial] = max_looker_resistance


            # this does not appear in the logic document, but it was in legacy code; not sure if we should retain?
            if DF_LOOKER[antimicrobial] != "R":
                for gene in drugs_to_genes[antimicrobial]:
                    if gene in self.LOW_DEPTH_OF_COVERAGE_LIST and gene not in self.GENES_WITH_VALID_DELETIONS:
                        DF_LOOKER[antimicrobial] = "Insufficient coverage in locus"
                        break


        # ENABLE FIELD RENAMING

        # as per rule 6.1.2 and 6.1.3, the lineage field is the main_lineage field from TBProfiler and the ID field is the same as ID in the LIMS report
        DF_LOOKER["lineage"] = LINEAGE
        DF_LOOKER["ID"] = LINEAGE_ENGLISH

        DF_LOOKER["analysis_date"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        DF_LOOKER["operator"] = OPERATOR

        DF_LOOKER.to_csv("{}.looker_report.csv".format(self.output_prefix), index=False)
