
import pandas as pd
import datetime

class Looker:
    """Class to create the Looker report based on the Laboratorian report and other inputs.
    
    Attributes:
        RESISTANCE_RANKING (dict[str, int]): A dictionary ranking the resistance annotations from 
            highest to lowest priority.
        logger (_type_): the logging instance
        OUTPUT_PREFIX (str): the output file prefix
        DF_LABORATORIAN (pd.DataFrame): the laboratorian report dataframe
        LOW_DEPTH_OF_COVERAGE_LIST (list[str]): list of genes with low depth of coverage
        GENES_WITH_VALID_DELETIONS (dict[str, list[int]]): dictionary of genes (to the deletion genomic positions) with valid deletions
        GENE_TO_ANTIMICROBIAL_DRUG_NAME (dict[str, dict[str, dict[str, str]]]): a
            dictionary mapping drugs to lims report columns to associated drugs and 
            their respective column names
    
    Methods:
        create_looker_report(SAMPLE_NAME: str, SEQUENCING_METHOD: str, LINEAGE: str, 
            LINEAGE_ENGLISH: str, OPERATOR: str) -> None: generates the looker report
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
    
    def __init__(self, logger, OUTPUT_PREFIX: str, DF_LABORATORIAN: pd.DataFrame, LOW_DEPTH_OF_COVERAGE_LIST: list[str], GENES_WITH_VALID_DELETIONS: dict[str, list[int]], GENE_TO_ANTIMICROBIAL_DRUG_NAME: dict[str, dict[str, dict[str, str]]]) -> None:
        """Initializes the Looker class.
        
        Args:
            logger (_type_): the logging instance
            OUTPUT_PREFIX (str): the output file prefix
            DF_LABORATORIAN (pd.DataFrame): the laboratorian report dataframe
            LOW_DEPTH_OF_COVERAGE_LIST (list[str]): list of genes with low depth of coverage
            GENES_WITH_VALID_DELETIONS (dict[str, list[int]]): dictionary of genes (to the deletion genomic positions) with valid deletions
            GENE_TO_ANTIMICROBIAL_DRUG_NAME (dict[str, dict[str, dict[str, str]]]): a
                dictionary mapping drugs to lims report columns to associated drugs and 
                their respective column names
        """    
        self.logger = logger
        self.OUTPUT_PREFIX = OUTPUT_PREFIX
        self.DF_LABORATORIAN = DF_LABORATORIAN
        self.LOW_DEPTH_OF_COVERAGE_LIST = LOW_DEPTH_OF_COVERAGE_LIST
        self.GENES_WITH_VALID_DELETIONS = GENES_WITH_VALID_DELETIONS
        self.GENE_TO_ANTIMICROBIAL_DRUG_NAME = GENE_TO_ANTIMICROBIAL_DRUG_NAME

    def create_looker_report(self, SAMPLE_NAME: str, SEQUENCING_METHOD: str, LINEAGE: str, LINEAGE_ENGLISH: str, OPERATOR: str) -> None:
        """Creates the Looker report; lists each antimicrobial and indicates the highest resistance level found

        Args:
            SAMPLE_NAME (str): the name of the sample
            SEQUENCING_METHOD (str): the method used to sequence the sample
            LINEAGE (str): the lineage as determined by TBProfiler
            LINEAGE_ENGLISH (str): the lineage rewritten to a human-readable format
            OPERATOR (str): the name of the person running the script
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
                
        antimicrobial_resistances = {}
        # iterate through laboratorian dataframe to extract highest mutation
        for antimicrobial in drugs_to_genes.keys():
            potential_looker_resistances = self.DF_LABORATORIAN[self.DF_LABORATORIAN["antimicrobial"] == antimicrobial]["looker_interpretation"]

            try:
                max_looker_resistance = max(potential_looker_resistances, key=lambda x: self.RESISTANCE_RANKING[x])
            except:
                max_looker_resistance = "NA"
                
            antimicrobial_resistances[antimicrobial] = max_looker_resistance

            # this does not appear in the logic document, but it was in legacy code; not sure if we should retain?
            if max_looker_resistance != "R":
                for gene in drugs_to_genes[antimicrobial]:
                    if gene in self.LOW_DEPTH_OF_COVERAGE_LIST and gene not in self.GENES_WITH_VALID_DELETIONS.keys():
                        antimicrobial_resistances[antimicrobial] = "Insufficient coverage in locus"
                        break

        for antimicrobial in sorted(antimicrobial_resistances.keys()):
            # order alphabetically; probably doesn't matter but makes it easy to compare to old versions of the script
            DF_LOOKER[antimicrobial] = antimicrobial_resistances[antimicrobial]
        
        # as per rule 6.1.2 and 6.1.3, the lineage field is the main_lineage field from TBProfiler and the ID field is the same as ID in the LIMS report
        DF_LOOKER["lineage"] = LINEAGE
        DF_LOOKER["ID"] = LINEAGE_ENGLISH

        DF_LOOKER["analysis_date"] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        DF_LOOKER["operator"] = OPERATOR

        DF_LOOKER.to_csv("{}.looker_report.csv".format(self.OUTPUT_PREFIX), index=False)
