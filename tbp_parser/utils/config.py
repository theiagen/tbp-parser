from typing import Any
import argparse
import logging
import yaml

logger = logging.getLogger(__name__)

class Configuration:
    """Configuration class for TBP Parser.

    This class handles all configuration settings and input parameters.

    Attributes:
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
    """
    # Shared type tracking across all Configuration instances
    _CONFIGURABLE_INPUTS = {}


    def __setattr__(self, name: str, value: Any) -> None:
        """Custom setattr to automatically track configurable inputs that are uppercase
        without impacting other attributes.

        Args:
            name (str): the name of the attribute to set
            value (any): the value to set the attribute to
        """
        if not name.startswith("_"):
            if name not in vars(self):
                # log only when setting new attributes or for the first time
                logger.debug(f"'{name}': {value}")

            if name.isupper():
                # track initial configurable inputs and their types enforced by argparse
                if name not in self._CONFIGURABLE_INPUTS:
                    self._CONFIGURABLE_INPUTS[name] = type(value)
                else:
                    # enforce type checking for subsequent assignments to configurable inputs (from config file)
                    expected_type = self._CONFIGURABLE_INPUTS[name]

                    # allow int/float interchangeability
                    if expected_type in (int, float) and isinstance(value, (int, float)):
                        value = expected_type(value)

                    if type(value) != expected_type:
                        raise TypeError(f"Type mismatch for attribute '{name}': expected {expected_type}, got {type(value)}")
        super().__setattr__(name, value)


    def __init__(self, options: argparse.Namespace) -> None:
        """Initialize the Configuration class

        Args:
            options (argparse.Namespace): an object with the input arguments provided at runtime
        """
        # INITIALIZE FIXED INPUTS (lowercase only)
        # main input files
        self.input_json = options.input_json
        self.input_bam = options.input_bam
        self.config = options.config
        # files to be parsed once and used across multiple classes
        self.tbdb_bed = options.tbdb_bed
        self.err_bed = options.err_bed
        self.gene_tier_tsv = options.gene_tier_tsv
        self.promoter_regions_tsv = options.promoter_regions_tsv
        self.lims_report_format_yml = options.lims_report_format_yml

        # INITIALIZE CONFIGURABLE INPUTS (always uppercase)
        # qc options
        self.MIN_DEPTH = options.min_depth
        self.MIN_PERCENT_COVERAGE = options.min_percent_coverage
        self.MIN_READ_SUPPORT = options.min_read_support
        self.MIN_FREQUENCY = options.min_frequency
        self.MIN_PERCENT_LOCI_COVERED = options.min_percent_loci_covered
        self.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY = options.do_not_treat_r_mutations_differently
        # text inputs
        self.SEQUENCING_METHOD = options.sequencing_method
        self.OUTPUT_PREFIX = options.output_prefix
        self.OPERATOR = options.operator
        self.FIND_AND_REPLACE = options.find_and_replace
        # tngs-specific options
        self.TNGS = options.tngs
        self.USE_ERR_AS_BRR = options.use_err_as_brr
        self.TNGS_READ_SUPPORT_BOUNDARIES = [int(x) for x in options.tngs_read_support_boundaries.split(",")]
        self.TNGS_FREQUENCY_BOUNDARIES = [float(x) for x in options.tngs_frequency_boundaries.split(",")]
        # logging options
        self.DEBUG = options.debug
        # tngs-specific qc options that are hold-overs from prior versions; retained for backwards compatibility
        self.TNGS_SPECIFIC_QC_OPTIONS = {
            "RRS_FREQUENCY": options.rrs_frequency,
            "RRS_READ_SUPPORT": options.rrs_read_support,
            "RRL_FREQUENCY": options.rrl_frequency,
            "RRL_READ_SUPPORT": options.rrl_read_support,
            "ETHA237_FREQUENCY": options.etha237_frequency,
            "RPOB449_FREQUENCY": options.rpob449_frequency,
        }

        # configuration file overwrite
        if self.config != "":
            logger.info("Overwriting variables with the provided config file")
            self.overwrite_variables()

        # NOT INPUTS/CONFIGURABLES
        # but will be populated later on and used across multiple classes
        # self._GENE_TO_ANTIMICROBIAL_DRUG_NAME = {}

    def overwrite_variables(self) -> None:
        """This function overwrites the input variables provided at runtime with those
        from the config file
        """
        with open(self.config, "r") as config:
            settings = yaml.safe_load(config)

            # only overwrite global variables if they are valid uppercase/overwritable attributes of self
            for key, value in settings.items():
                if key in self._CONFIGURABLE_INPUTS:
                    setattr(self, key, value)
                    logger.debug(f"'{key}': {value}")