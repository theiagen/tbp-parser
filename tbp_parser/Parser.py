import logging
import subprocess
import globals as globals_
import sys
import re
import yaml
from Coverage import Coverage
from Laboratorian import Laboratorian
from Looker import Looker
from LIMS import LIMS


class Parser:
    """This class orchestrates the different modules within the tbp_parser tool.
    
    It initializes with the input arguments provided at runtime, creates the standard
    look-up dictionaries, checks for dependencies, and runs the main logic to generate
    the desired reports.
    """

    def __init__(self, options):
        """Initialize the Parser class

        Args:
            options (argparse.NameSpace): an object with the input arguments provided at runtime
        """        
        logging.basicConfig(encoding='utf-8', level=logging.ERROR, stream=sys.stderr)
        self.logger = logging.getLogger(__name__)
        self.input_json = options.input_json
        self.input_bam = options.input_bam
        self.config = options.config
        self.verbose = options.verbose
        self.debug = options.debug
        self.output_prefix = options.output_prefix
        
        self.MIN_PERCENT_COVERAGE = options.min_percent_coverage
        self.MIN_DEPTH = options.min_depth
        
        # files used to create the standard dictionaries
        self.tbdb_bed = options.tbdb_bed
        self.gene_tier_file = options.gene_tier_file
        self.promoter_regions_file = options.promoter_regions_file
        
        # reevaluate the following global variables -- they don't need to be globals
        self.TNGS = options.tngs
        self.SEQUENCING_METHOD = options.sequencing_method
        self.MIN_READ_SUPPORT = options.min_read_support
        self.MIN_FREQUENCY = options.min_frequency
        self.MIN_LOCUS_PERCENTAGE = options.min_percent_locus_covered
        """The minimum percentage of LIMS genes to pass QC for MTBC identification to occur"""
        
        self.OPERATOR = options.operator

        # this could cause issues if someone does more than one comma, but in that case, they deserve the error
        globals_.TNGS_READ_SUPPORT_BOUNDARIES = [int(x) for x in options.tngs_read_support_boundaries.split(",")]
        globals_.TNGS_FREQUENCY_BOUNDARIES = [float(x) for x in options.tngs_frequency_boundaries.split(",")]
        globals_.TREAT_R_AS_S = options.treat_r_mutations_as_s

        if self.verbose:
            self.logger.setLevel(logging.INFO)
            self.logger.info("PARSER:__init__:Verbose mode enabled")
        else:
            self.logger.setLevel(logging.ERROR)

        if self.debug:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug("PARSER:__init__:Debug mode enabled")

        if self.config != "":
            self.logger.info("PARSER:__init__:Overwriting variables with the provided config file")
            self.overwrite_variables()

    def create_standard_dictionaries(self) -> tuple[dict, dict, dict, dict, dict]:
        """Parses the provided coverage regions bed file to create the standard
        look-up dictionaries used throughout tbp-parser
        
        Also, if tNGS is enabled, creates the tNGS regions dictionary to
        confirm that the mutations are within the expected regions [tNGS only]

        For genes with suffixes like gene_1, gene_2, they are consolidated under the
        common gene name, with each entry stored as a sub-item in a nested dictionary.
        
        That is, for a gene named "gene" with split regions "gene_1" and "gene_2", and
        a gene named "gene3" without splits, the structure will be:
        
        globals_.TNGS_REGIONS = {
            "gene": {
                "gene_1": [start_pos_1, end_pos_1],
                "gene_2": [start_pos_2, end_pos_2]
            },
            "gene3": [start_pos_3, end_pos_3]
            ...
        }

        Returns:
            tuple[dict, dict, dict, dict, dict]: gene to antimicrobial drug name dictionary,
            gene to locus tag dictionary, tNGS regions dictionary, gene to tier dictionary,
            and promoter regions dictionary
        """        
        GENE_TO_ANTIMICROBIAL_DRUG_NAME = {}
        GENE_TO_LOCUS_TAG = {}
        TNGS_REGIONS = {}
        with open(self.tbdb_bed, 'r') as bed_file:
            for line in bed_file:
                cols = line.strip().split('\t')
                gene_name = cols[4]
                drugs = cols[5]
                locus_tag = cols[3]
                                
                GENE_TO_ANTIMICROBIAL_DRUG_NAME[gene_name] = drugs
                GENE_TO_LOCUS_TAG[gene_name] = locus_tag

                if self.tngs:
                    start_pos = cols[1]
                    end_pos = cols[2]
                    
                    # check if primer is split
                    match = re.match(r'(.+)_(\d+)$', gene_name)
                    if match:
                        base_gene_name = match.group(1)
                        if base_gene_name not in TNGS_REGIONS:
                            TNGS_REGIONS[base_gene_name] = {}
                        TNGS_REGIONS[base_gene_name][gene_name] = [start_pos, end_pos]

                    else:
                        TNGS_REGIONS[gene_name] = [start_pos, end_pos]

        GENE_TO_TIER = {}
        with open(self.gene_tier_file, 'r') as tier_file:
            for line in tier_file:
                cols = line.strip().split('\t')
                gene_name = cols[0]
                tier = cols[1]
                
                GENE_TO_TIER[gene_name] = tier
       
        PROMOTER_REGIONS = {}
        with open(self.promoter_regions_file, 'r') as promoter_file:
            for line in promoter_file:
                cols = line.strip().split('\t')
                gene_name = cols[0]
                start_pos = cols[1]
                end_pos = cols[2]
                
                # check if gene has a second promoter in cols[3] and cols[4]
                try:
                    if cols[3] != "" and cols[4] != "":
                        start_pos_2 = cols[3]
                        end_pos_2 = cols[4]
                        PROMOTER_REGIONS[gene_name] = [[start_pos, end_pos], [start_pos_2, end_pos_2]]
                        continue
                except IndexError:
                    # no second promoter
                    pass
                
                # avoid overwriting if already exists (if a gene has multiple promoters)
                if gene_name not in PROMOTER_REGIONS:
                    PROMOTER_REGIONS[gene_name] = [start_pos, end_pos]
                                
        return GENE_TO_ANTIMICROBIAL_DRUG_NAME, GENE_TO_LOCUS_TAG, TNGS_REGIONS, GENE_TO_TIER, PROMOTER_REGIONS

    def overwrite_variables(self) -> None:
        """This function overwrites the input variables provided at runtime with those 
        from the config file
        """
        with open(self.config, "r") as config:
            settings = yaml.safe_load(config)

            for key, value in settings.items():
                if key.replace("self.", "") in vars(self):
                    setattr(self, key.replace("self.", ""), value)
                    self.logger.info("PARSER:overwrite_variables:self.{} has been overwritten with a config-specified value".format(key))
                if key.replace("globals.", "") in dir(globals_):
                    setattr(globals_, key.replace("globals.", ""), value) 
                    self.logger.info("PARSER:overwrite_variables:globals.{} has been overwritten with a config-specified value".format(key))

    def check_dependency_exists(self) -> None:
        """This function confirms that samtools is installed and available; if it is
        not, the program exits with an error message.
        """
        result = subprocess.run(
            ["samtools", "--version"],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE,
          text=True,
        )

        if result.returncode != 0:
            self.logger.critical("PARSER:Error: samtools not found. Please install samtools and try again.")
            sys.exit(1)

    def run(self) -> None:
        """This function runs the main logic of the Parser class, orchestrating the
        different modules to generate the desired reports.
        """    
        self.logger.info("PARSER:run:Checking for dependencies")
        self.check_dependency_exists()

        GENE_TO_ANTIMICROBIAL_DRUG_NAME, GENE_TO_LOCUS_TAG, TNGS_REGIONS, GENE_TO_TIER, PROMOTER_REGIONS = self.create_standard_dictionaries()
        
        self.logger.info("PARSER:run:Calculating coverage statistics")
        coverage = Coverage(self.logger, self.input_bam, self.output_prefix, self.coverage_regions)
        COVERAGE_DICTIONARY, AVERAGE_LOCI_COVERAGE, LOW_DEPTH_OF_COVERAGE_LIST = coverage.get_coverage(self.MIN_PERCENT_COVERAGE, self.MIN_DEPTH)

        # theoretically code should work up until here

        self.logger.info("PARSER:run:Creating Laboratorian report")
        laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix, 
                                    self.MIN_DEPTH, self.MIN_FREQUENCY, self.MIN_READ_SUPPORT, 
                                    COVERAGE_DICTIONARY, LOW_DEPTH_OF_COVERAGE_LIST, GENE_TO_ANTIMICROBIAL_DRUG_NAME, GENE_TO_LOCUS_TAG, TNGS_REGIONS, GENE_TO_TIER, PROMOTER_REGIONS, self.TNGS)
        DF_LABORATORIAN = laboratorian.create_laboratorian_report()

        self.logger.info("PARSER:run:Creating LIMS report")
        lims = LIMS(self.logger, self.input_json, self.output_prefix, LOW_DEPTH_OF_COVERAGE_LIST, laboratorian.SAMPLE_NAME, DF_LABORATORIAN, laboratorian.positional_qc_fails, laboratorian.genes_with_valid_deletions)
        lims.create_lims_report(self.TNGS, self.MIN_LOCUS_PERCENTAGE, self.OPERATOR)

        self.logger.info("PARSER:run:Creating Looker report")
        looker = Looker(self.logger, self.output_prefix, DF_LABORATORIAN, LOW_DEPTH_OF_COVERAGE_LIST, laboratorian.genes_with_valid_deletions, GENE_TO_ANTIMICROBIAL_DRUG_NAME)
        looker.create_looker_report(laboratorian.SAMPLE_NAME, self.SEQUENCING_METHOD, lims.LINEAGE, lims.LINEAGE_ENGLISH, self.OPERATOR)

        self.logger.info("PARSER:run:Creating coverage report")
        coverage.create_coverage_report(COVERAGE_DICTIONARY, AVERAGE_LOCI_COVERAGE, laboratorian.genes_with_valid_deletions)

        # DO MASSIVE RENAMING HERE!!!!!

        self.logger.info("PARSER:run:Parsing completed")