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
    """This class orchestrates the different modules within the tbp_parser tool."""

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
        self.coverage_regions = options.coverage_regions
        
        # reevaluate the following global variables
        globals_.TNGS = options.tngs
        globals_.TREAT_R_AS_S = options.treat_r_mutations_as_s
        globals_.MIN_DEPTH = options.min_depth
        globals_.COVERAGE_THRESHOLD = options.min_percent_coverage
        globals_.SEQUENCING_METHOD = options.sequencing_method
        globals_.MIN_READ_SUPPORT = options.min_read_support
        globals_.MIN_FREQUENCY = options.min_frequency
        globals_.MIN_LOCUS_PERCENTAGE = options.min_percent_locus_covered
        
        # TO-DO: consider dropping these parameters?
        globals_.RRS_FREQUENCY = options.rrs_frequency
        globals_.RRS_READ_SUPPORT = options.rrs_read_support
        globals_.RRL_FREQUENCY = options.rrl_frequency
        globals_.RRL_READ_SUPPORT = options.rrl_read_support
        globals_.RPOB449_FREQUENCY = options.rpob449_frequency
        globals_.ETHA237_FREQUENCY = options.etha237_frequency
        ###
        globals_.OPERATOR = options.operator

        # this could cause issues if someone does more than one comma, but in that case, they deserve the error
        globals_.TNGS_READ_SUPPORT_BOUNDARIES = [int(x) for x in options.tngs_read_support_boundaries.split(",")]
        globals_.TNGS_FREQUENCY_BOUNDARIES = [float(x) for x in options.tngs_frequency_boundaries.split(",")]

        if self.verbose:
            self.logger.setLevel(logging.INFO)
            self.logger.info("PARSER:__init__:Verbose mode enabled")
        else:
            self.logger.setLevel(logging.ERROR)

        if self.debug:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug("PARSER:__init__:Debug mode enabled")

        if globals_.TNGS:
            self.logger.debug("PARSER:__init__:Setting the tNGS regions dictionary")
            self.convert_bed_into_dictionary()

        if self.config != "":
            self.logger.info("PARSER:__init__:Overwriting variables with the provided config file")
            self.overwrite_variables()

    def convert_bed_into_dictionary(self) -> None:
        """This function converts the `coverage_regions` bed file into a dictionary to
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
        """ 
        with open(self.coverage_regions, 'r') as bed_file:
            for line in bed_file:
                cols = line.strip().split('\t')

                start_pos = int(cols[1])
                end_pos = int(cols[2])
                gene_name = cols[4]

                # check if primer is split
                match = re.match(r'(.+)_(\d+)$', gene_name)
                if match:
                    base_gene_name = match.group(1)
                    if base_gene_name not in globals_.TNGS_REGIONS:
                        globals_.TNGS_REGIONS[base_gene_name] = {}
                    globals_.TNGS_REGIONS[base_gene_name][gene_name] = [start_pos, end_pos]

                else:
                    globals_.TNGS_REGIONS[gene_name] = [start_pos, end_pos]

            self.logger.debug("PARSER:convert_bed_into_dictionary:Finished processing coverage regions; {}".format(globals_.TNGS_REGIONS))

    def overwrite_variables(self) -> None:
        """This function overwrites the input variables provided at runtime with those from the config file"""
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
        """This function confirms that samtools is installed and available"""
        result = subprocess.run(
            ["samtools", "--version"],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE,
          text=True,
        )

        if result.returncode != 0:
            self.logger.critical("PARSER:Error: samtools not found. Please install samtools and try again.")
            sys.exit(1)

    def run(self):
        """
        This function runs the parsing module for the tb_parser tool.
        """    
        self.logger.info("PARSER:run:Checking for dependencies")
        self.check_dependency_exists

        self.logger.info("PARSER:run:Calculating coverage statistics")
        coverage = Coverage(self.logger, self.input_bam, self.output_prefix, self.coverage_regions)
        coverage.get_coverage()

        self.logger.info("PARSER:run:Creating Laboratorian report")
        laboratorian = Laboratorian(self.logger, self.input_json, self.output_prefix)
        laboratorian.create_laboratorian_report()

        self.logger.info("PARSER:run:Creating LIMS report")
        lims = LIMS(self.logger, self.input_json, self.output_prefix)
        lims.create_lims_report()

        self.logger.info("PARSER:run:Creating Looker report")
        looker = Looker(self.logger, self.output_prefix)
        looker.create_looker_report()

        self.logger.info("PARSER:run:Finalizing coverage report")
        coverage.create_coverage_report()

        self.logger.info("PARSER:run:Parsing completed")