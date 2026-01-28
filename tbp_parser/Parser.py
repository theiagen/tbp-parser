import logging
import subprocess
import globals as globals_
import sys
import fileinput
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
    
    Attributes:
        logger (logging.Logger): the logger object for logging messages
        input_json (str): the JSON file produced by TBProfiler
        input_bam (str): the BAM file produced by TBProfiler
        config (str): the configuration file to use, in YAML format        
        
        OPERATOR (str): the operator who ran the sequencing
        OUTPUT_PREFIX (str): the prefix for output files
        SEQUENCING_METHOD (str): the sequencing method used to generate the data
        
        tbdb_bed (str): the BED file containing the genes of interest, their locus tags, their associated antimicrobial, and their regions for QC calculations
        gene_tier_tsv (str): the TSV file mapping genes to their tier
        promoter_regions_tsv (str): the TSV file containing the promoter regions to include in interpretation designations
        
        MIN_PERCENT_COVERAGE (float): the minimum percentage of a region that has depth above the threshold set by MIN_DEPTH to pass QC 
        MIN_PERCENT_LOCI_COVERED (float): the minimum percentage of loci/genes in the LIMS report that must pass coverage QC for the sample to be identified as MTBC
        MIN_DEPTH (int): the minimum depth of coverage for a site to pass QC
        MIN_FREQUENCY (float): the minimum frequency for a mutation to pass QC
        MIN_READ_SUPPORT (int): the minimum read support for a mutation to pass QC
        
        TNGS (bool): whether tNGS mode is enabled
        DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY (bool): whether R mutations should be treated the same as S/U mutations for locus QC
        TNGS_READ_SUPPORT_BOUNDARIES (list[int]): the read support boundaries for tNGS QC reporting
        TNGS_FREQUENCY_BOUNDARIES (list[float]): the frequency boundaries for tNGS QC reporting
        err_bed (str | None): an optional BED file containing ranges that are essential for resistance [tNGS only]
        TNGS_SPECIFIC_QC_OPTIONS (dict): tNGS-specific QC options that are hold-overs from prior versions; retained for backwards compatibility
        USE_ERR_AS_BRR (bool): whether to use ERR regions in place of TBDB regions for breadth of coverage calculations [tNGS only]
        
    Methods:
        create_standard_dictionaries() -> tuple[dict[str, list[str]], dict[str, str], dict[str, list[int] | list[list[int]]], dict[str, str], dict[str, list[int] | list[list[int]]]]:
            Creates the standard look-up dictionaries used throughout tbp-parser
        overwrite_variables() -> None:
            Overwrites input variables with those from the config file
        check_dependency_exists() -> None:
            Checks if samtools is installed and available
        run() -> None:
            Runs the main logic of the Parser class to generate the reports
    
    """

    def __init__(self, options) -> None:
        """Initialize the Parser class

        Args:
            options (argparse.NameSpace): an object with the input arguments provided at runtime
        """        
        logging.basicConfig(encoding='utf-8', level=logging.ERROR, stream=sys.stderr)
        self.logger = logging.getLogger(__name__)
        self.input_json = options.input_json
        self.input_bam = options.input_bam
        
        self.config = options.config
        
        # test inputs
        self.SEQUENCING_METHOD = options.sequencing_method
        self.OUTPUT_PREFIX = options.output_prefix
        self.OPERATOR = options.operator
        
        # files used to create the standard dictionaries
        self.tbdb_bed = options.tbdb_bed
        self.gene_tier_tsv = options.gene_tier_tsv
        self.promoter_regions_tsv = options.promoter_regions_tsv
        
        self.MIN_PERCENT_COVERAGE = options.min_percent_coverage
        self.MIN_PERCENT_LOCI_COVERED = options.min_percent_loci_covered
        """The minimum percentage of LIMS genes to pass QC for MTBC identification to occur"""
        self.MIN_DEPTH = options.min_depth
        self.MIN_FREQUENCY = options.min_frequency
        self.MIN_READ_SUPPORT = options.min_read_support
        
        # tngs-specific options
        self.TNGS = options.tngs
        self.TNGS_READ_SUPPORT_BOUNDARIES = [int(x) for x in options.tngs_read_support_boundaries.split(",")]
        self.TNGS_FREQUENCY_BOUNDARIES = [float(x) for x in options.tngs_frequency_boundaries.split(",")]
        self.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY = options.do_not_treat_r_mutations_differently
        self.err_bed = options.err_bed
        self.USE_ERR_AS_BRR = options.use_err_as_brr

        # tngs-specific qc options that are hold-overs from prior versions; retained for backwards compatibility
        self.TNGS_SPECIFIC_QC_OPTIONS = {
            "RRS_FREQUENCY": options.rrs_frequency,
            "RRS_READ_SUPPORT": options.rrs_read_support,
            "RRL_FREQUENCY": options.rrl_frequency,
            "RRL_READ_SUPPORT": options.rrl_read_support,
            "ETHA237_FREQUENCY": options.etha237_frequency,
            "RPOB449_FREQUENCY": options.rpob449_frequency,
        }

        # logging options
        if options.verbose:
            self.logger.setLevel(logging.INFO)
            self.logger.info("PARSER:__init__:Verbose mode enabled")
        elif options.debug:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug("PARSER:__init__:Debug mode enabled")
        else:
            self.logger.setLevel(logging.ERROR)

        # configuration file overwrite
        if self.config != "":
            self.logger.info("PARSER:__init__:Overwriting variables with the provided config file")
            self.overwrite_variables()

    def create_standard_dictionaries(self) -> tuple[dict[str, list[str]], dict[str, str], dict[str, list[int] | dict[str, list[int]]], dict[str, str], dict[str, list[int] | list[list[int]]]]:
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
                             
                if self.TNGS:
                    start_pos = int(cols[1])
                    end_pos = int(cols[2])
                    
                    # check if primer is split
                    match = re.match(r'(.+)_(\d+)$', gene_name)
                    if match:
                        base_gene_name = match.group(1)
                        if base_gene_name not in TNGS_REGIONS:
                            TNGS_REGIONS[base_gene_name] = {}
                        TNGS_REGIONS[base_gene_name][gene_name] = [start_pos, end_pos]
                    else:
                        TNGS_REGIONS[gene_name] = [start_pos, end_pos]
                    
                    gene_name = base_gene_name if match else gene_name
                    
                GENE_TO_ANTIMICROBIAL_DRUG_NAME[gene_name] = drugs.split(",")
                GENE_TO_LOCUS_TAG[gene_name] = locus_tag

        GENE_TO_TIER = {}
        with open(self.gene_tier_tsv, 'r') as tier_file:
            for line in tier_file:
                cols = line.strip().split('\t')
                gene_name = cols[0]
                tier = cols[1]
                
                GENE_TO_TIER[gene_name] = tier
       
        PROMOTER_REGIONS = {}
        with open(self.promoter_regions_tsv, 'r') as promoter_file:
            for line in promoter_file:
                cols = line.strip().split('\t')
                gene_name = cols[0]
                start_pos = int(cols[1])
                end_pos = int(cols[2])
                
                # check if gene has a second promoter in cols[3] and cols[4]
                try:
                    if cols[3] != "" and cols[4] != "":
                        start_pos_2 = int(cols[3])
                        end_pos_2 = int(cols[4])
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
        coverage = Coverage(self.logger, self.input_bam, self.OUTPUT_PREFIX, self.tbdb_bed, TNGS_REGIONS, self.USE_ERR_AS_BRR, self.err_bed)
        COVERAGE_DICTIONARY, LOW_DEPTH_OF_COVERAGE_LIST = coverage.get_coverage(self.MIN_PERCENT_COVERAGE, self.MIN_DEPTH, self.TNGS)

        self.logger.info("PARSER:run:Creating Laboratorian report")
        laboratorian = Laboratorian(self.logger, self.input_json, self.OUTPUT_PREFIX, 
                                    self.MIN_DEPTH, self.MIN_FREQUENCY, self.MIN_READ_SUPPORT, 
                                    COVERAGE_DICTIONARY, LOW_DEPTH_OF_COVERAGE_LIST, GENE_TO_ANTIMICROBIAL_DRUG_NAME, 
                                    GENE_TO_LOCUS_TAG, TNGS_REGIONS, GENE_TO_TIER, PROMOTER_REGIONS, self.TNGS, 
                                    self.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY, self.TNGS_READ_SUPPORT_BOUNDARIES,
                                    self.TNGS_FREQUENCY_BOUNDARIES, self.TNGS_SPECIFIC_QC_OPTIONS)
        DF_LABORATORIAN = laboratorian.create_laboratorian_report()

        self.logger.info("PARSER:run:Creating LIMS report")
        lims = LIMS(self.logger, self.input_json, self.OUTPUT_PREFIX, LOW_DEPTH_OF_COVERAGE_LIST, 
                    laboratorian.SAMPLE_NAME, DF_LABORATORIAN, laboratorian.positional_qc_fails,
                    laboratorian.genes_with_valid_deletions)
        lims.create_lims_report(self.TNGS, self.MIN_PERCENT_LOCI_COVERED, self.OPERATOR)

        self.logger.info("PARSER:run:Creating Looker report")
        looker = Looker(self.logger, self.OUTPUT_PREFIX, DF_LABORATORIAN, LOW_DEPTH_OF_COVERAGE_LIST, laboratorian.genes_with_valid_deletions, GENE_TO_ANTIMICROBIAL_DRUG_NAME)
        looker.create_looker_report(laboratorian.SAMPLE_NAME, self.SEQUENCING_METHOD, lims.LINEAGE, lims.LINEAGE_ENGLISH, self.OPERATOR)

        self.logger.info("PARSER:run:Creating coverage report")
        coverage.create_coverage_report(laboratorian.SAMPLE_NAME, laboratorian.genes_with_valid_deletions)

        if len(globals_.OUTPUT_RENAMING) > 0:
            self.logger.info("PARSER:run:Renaming output columns as specified in globals.OUTPUT_RENAMING")
            
            # only rename the fields in the final output files
            file_suffixes = [".lims_report.csv", ".laboratorian_report.csv", ".looker_report.csv", ".coverage_report.csv"]
            for output_suffix in file_suffixes:
                for line in fileinput.input(self.OUTPUT_PREFIX + output_suffix, inplace=True):
                    # match full words to prevent partial replacements
                    print(re.sub('|'.join(r'\b%s\b' % re.escape(original) for original in globals_.OUTPUT_RENAMING.keys()),
                                 lambda match: globals_.OUTPUT_RENAMING[match.group(0)], line), end="")
                    
        self.logger.info("PARSER:run:Parsing completed")
