
import logging

from typing import Dict
from utils.config import Configuration
from utils.helper import Helper
from variant import Variant
from coverage import LocusCoverage

logger = logging.getLogger(__name__)

class VariantQC:

    def __init__(self, config: Configuration, ):
        self.config = config

    @staticmethod
    def is_deletion(variant: Variant) -> bool:
        return "del" in variant.nucleotide_change

    def fails_tngs_specific_qc(self, variant: Variant) -> bool:
        """Checks if a mutation (tNGS only) fails the tNGS-specific QC checks

        Args:
            TNGS_SPECIFIC_QC_OPTIONS (dict[str, float]): a dictionary containing tNGS specific QC options

        Returns:
            bool: true if the mutation fails QC, false if the mutation passes QC
        """
        RRS_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRS_FREQUENCY"]
        RRS_READ_SUPPORT = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRS_READ_SUPPORT"]
        RRL_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRL_FREQUENCY"]
        RRL_READ_SUPPORT = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRL_READ_SUPPORT"]
        ETHA237_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["ETHA237_FREQUENCY"]
        RPOB449_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["RPOB449_FREQUENCY"]

        if variant.gene_name == "rrs":
            if variant.freq < RRS_FREQUENCY or variant.read_support < RRS_READ_SUPPORT:
                return True
        elif variant.gene_name == "rrl":
            if variant.freq < RRL_FREQUENCY or variant.read_support < RRL_READ_SUPPORT:
                return True
        elif variant.gene_name == "ethA" and Helper.get_position(variant.protein_change) == 237:
            if variant.freq < ETHA237_FREQUENCY:
                return True
        elif variant.gene_name == "rpoB" and Helper.get_position(variant.protein_change) == 449:
            if variant.freq < RPOB449_FREQUENCY:
                return True
        return False

    def fails_tngs_boundary_qc(self, variant: Variant) -> bool:
        """Checks if a mutation (tNGS only) fails the boundary QC checks

        Args:
            TNGS_READ_SUPPORT_BOUNDARIES (list[int]): the lower and upper read support boundaries
            TNGS_FREQUENCY_BOUNDARIES (list[float]): the lower and upper frequency boundaries

        Returns:
            bool: true if the mutation fails QC, false if the mutation passes QC
        """
        lower_rs = self.config.TNGS_READ_SUPPORT_BOUNDARIES[0]
        upper_rs = self.config.TNGS_READ_SUPPORT_BOUNDARIES[1]
        lower_f = self.config.TNGS_FREQUENCY_BOUNDARIES[0]
        upper_f = self.config.TNGS_FREQUENCY_BOUNDARIES[1]

        FAILS_QC = False

        if (lower_rs <= variant.read_support and variant.read_support < upper_rs):
            if (variant.freq < upper_f):
                FAILS_QC = True
        elif (variant.read_support >= upper_rs):
            if (variant.freq < lower_f):
                FAILS_QC = True
        elif (variant.read_support < lower_rs):
            FAILS_QC = True

        if FAILS_QC:
            logger.debug(f"{variant} at [RS:{variant.read_support}, F:{variant.freq}] fails tNGS boundary QC: [RS:({lower_rs},{upper_rs}), F:({lower_f},{upper_f})]")

        return FAILS_QC

    def add_qc_warning(self, variants: list[Variant], locus_coverage_map: Dict[str, LocusCoverage]) -> None:
        """Adds QC warnings if a mutation either has poor positional quality or locus quality
          Assigns `fails_qc` and `warning` attributes to each Variant object if applicable.
        
          Args:
            variants (list[Variant]): List of Variant objects to check for QC warnings

          Returns:
            None
        """
        for variant in variants:
            locus_coverage = locus_coverage_map.get(variant.gene_id, None)
            if (
                variant.depth < self.config.MIN_DEPTH or
                variant.freq < self.config.MIN_FREQUENCY or
                variant.read_support < self.config.MIN_READ_SUPPORT or
                self.config.TNGS and self.fails_tngs_specific_qc(variant)
            ):
                variant.fails_qc = True
                variant.warning.add("Failed quality in the mutation position")

            elif VariantQC.is_deletion(variant):
                if (
                    variant.depth > 0 and
                    variant.depth < self.config.MIN_DEPTH
                ):
                    # 4.2.1.2 - postiional qc fail, deletion with some depth but not enough
                    variant.fails_qc = True
                    variant.warning.add("Failed quality in the mutation position")

                elif (
                    variant.depth == 0 and
                    variant.freq >= self.config.MIN_FREQUENCY
                ):
                    # 4.2.1.3 - deletion with zero depth but good frequency
                    pass

                elif (variant.freq < self.config.MIN_FREQUENCY):
                    # frequency is poor -- positional qc fail
                    variant.fails_qc = True
                    variant.warning.add("Failed quality in the mutation position")

            if self.config.TNGS:
                if self.fails_tngs_boundary_qc(variant):
                    variant.fails_qc = True
                    variant.warning.add("Failed quality in the mutation position")

                # is mutation outside of tNGS primer regions?
                if (
                  not locus_coverage or
                  not locus_coverage.contains_position(variant.pos)
                ):
                    variant.fails_qc = True
                    variant.warning.add("This mutation is outside the expected region")
                    # i feel like this should be moved somewhere else...
                    variant.rationale = "NA"
                    variant.confidence = "NA"
                    variant.looker_interpretation = "NA"
                    variant.mdl_interpretation = "NA"

            # checking locus qc now
            if (
                locus_coverage and
                locus_coverage.has_breadth_below(self.config.MIN_PERCENT_COVERAGE)
            ):
                if (
                    VariantQC.is_deletion(variant) and # 4.2.2.2 - locus qc fail and a deletion
                    variant.fails_qc
                ):
                    # this mutation also failed positional qc, so we need to add the locus warning
                    variant.warning.add("Insufficient coverage in locus")

                else: # 4.2.2.3 - locus qc fail but not a deletion
                    if variant.mdl_interpretation == "R" and not self.config.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY:
                        if not variant.fails_qc:
                            # 4.2.2.3.3 - R mutation with locus qc fail but NOT positional qc fail; add warning DO NOT not overwrite interpretation
                            variant.warning.add("Insufficient coverage in locus")
                        else:
                            # 4.2.2.3.4 - R mutation with BOTH positional and locus qc fail; add warning and overwrite interpretation
                            variant.warning.add("Insufficient coverage in locus")
                            variant.looker_interpretation = "Insufficient Coverage"
                            variant.mdl_interpretation = "Insufficient Coverage"

                    elif variant.mdl_interpretation == "U" or variant.mdl_interpretation == "S":
                        # 4.2.2.3.2 - non-R mutation with locus qc fail; add warning and overwrite interpretation (or if indicated that r mutations should be treated the same as s/u)
                        variant.warning.add("Insufficient coverage in locus")
                        variant.looker_interpretation = "Insufficient Coverage"
                        variant.mdl_interpretation = "Insufficient Coverage"

            if (
                not variant.fails_qc and
                VariantQC.is_deletion(variant)
            ):
                variant.is_valid_deletion = True

    @staticmethod
    def get_genes_with_valid_deletions(variants: list[Variant]) -> dict:
        """Returns a dict of gene names to lists of genomic positions for deletions
        that passed QC checks.

        Args:
            variants: List of Variant objects (after QC warnings have been applied)

        Returns:
            dict[str, list[int]]: gene_name -> list of genomic positions with valid deletions
        """
        genes_with_valid_deletions = {}
        for variant in variants:
            if getattr(variant, 'is_valid_deletion', False):
                positions = Helper.get_mutation_genomic_positions(variant.pos, variant.nucleotide_change)
                if variant.gene_name not in genes_with_valid_deletions:
                    genes_with_valid_deletions[variant.gene_name] = []
                genes_with_valid_deletions[variant.gene_name].extend(positions if isinstance(positions, list) else [positions])
        return genes_with_valid_deletions
