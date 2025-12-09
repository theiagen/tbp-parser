
import globals as globals_
import pandas as pd
import datetime

class Looker:
    """ 
    This class creates the CDPH Looker report.

    It has one function:
        - create_looker_report: creates the Looker report CSV file 
    """
    def __init__(self, logger, output_prefix):
        self.logger = logger
        self.output_prefix = output_prefix

    def create_looker_report(self):
        """
        This function recieves the input json and laboratorian report to
        write the Looker report that includes the following information: 
            - sample_id: the sample name
            - for each antimicrobial, indication if resistant (R) or susceptible (S)
        """
        self.logger.info("LOOKER:Within Looker class create_looker_report function")

        DF_LOOKER = pd.DataFrame({
            "sample_id": globals_.SAMPLE_NAME, 
          "output_seq_method_type": globals_.SEQUENCING_METHOD
            }, index=[0])

        # A dictionary that ranks the resistances on severity
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
        
        # iterate through laboratorian dataframe to extract highest mutation
        for antimicrobial in globals_.ANTIMICROBIAL_DRUG_NAME_LIST:
            self.logger.debug("LOOKER:Now extracting the highest mutation ranking for this antimicrobial: {}".format(antimicrobial))
            potential_looker_resistances = globals_.DF_LABORATORIAN[globals_.DF_LABORATORIAN["antimicrobial"] == antimicrobial]["looker_interpretation"]

            # this is a crazy one liner:
            # basically, it gets the max resistance ranking (R > R-Interim > U > S-Interim > S) for all resistance annotations for a drug
            try:
                max_looker_resistance = [annotation for annotation, rank in RESISTANCE_RANKING.items() if rank == max([RESISTANCE_RANKING[interpretation] for interpretation in potential_looker_resistances])]
            except:
                max_looker_resistance = ["NA"]
            DF_LOOKER[antimicrobial] = max_looker_resistance[0]
            self.logger.debug("LOOKER:The max Looker resistance for this antimicrobial is {}".format(max_looker_resistance[0]))

            for gene in globals_.ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME[antimicrobial]:
                # indicate warning if any genes failed to achieve 100% coverage_threshold
                if DF_LOOKER[antimicrobial][0] != "R" and gene in globals_.LOW_DEPTH_OF_COVERAGE_LIST and gene not in globals_.GENES_WITH_DELETIONS:
                    DF_LOOKER[antimicrobial] = "Insufficient coverage in locus"

        # get time stamp
        current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')

        # as per rule 6.1.2 and 6.1.3, the lineage field is the main_lineage field from TBProfiler and the ID field is the same as ID in the LIMS report
        DF_LOOKER["lineage"] = globals_.LINEAGE
        DF_LOOKER["ID"] = globals_.LINEAGE_ENGLISH

        DF_LOOKER["analysis_date"] = current_time
        DF_LOOKER["operator"] = globals_.OPERATOR

        # write to file
        DF_LOOKER.to_csv("{}.looker_report.csv".format(self.output_prefix), index=False)

        self.logger.info("LOOKER:Looker report created, now exiting function\n")