
import logging
from pydantic import BaseModel, Field
from typing import Dict, Optional
from utils import Configuration, Helper
from variant import Variant
from coverage import LocusCoverage, TargetCoverage
logger = logging.getLogger(__name__)
class QCResult(BaseModel):
    """Result of a QC check on a variant. Includes Variant attributes to update/overwrite based on the QC outcome."""
    type: Optional[str] = None
    nucleotide_change: Optional[str] = None
    protein_change: Optional[str] = None
    rationale: Optional[str] = None
    confidence: Optional[str] = None
    looker_interpretation: Optional[str] = None
    mdl_interpretation: Optional[str] = None
    fails_qc: bool = False
    warning: set[str] = Field(default_factory=set)

    def repr_filtered(self, attr_list: list[str]) -> str:
        """Return a string representation of the QCResult, filtered to only include specified attributes."""
        attrs = {k: v for k, v in self.model_dump(exclude_none=True).items() if k in attr_list}
        return f"QCResult({', '.join([f'{k}={v}' for k, v in attrs.items()])})"

class VariantQC:

    POSITIONAL_QC_WARNING = "Failed quality in the mutation position"
    LOCUS_QC_WARNING = "Insufficient coverage in locus"

    def __init__(self, config: Configuration):
        self.config = config

    def apply_qc(
        self,
        variants: list[Variant],
        locus_coverage_map: Dict[str, LocusCoverage],
        target_coverage_map: Dict[str, TargetCoverage],
        genes_with_valid_deletions: dict[str, list[Variant]],
    ) -> list[Variant]:
        """Apply all QC checks to a list of variants.

        Iterates through variants and applies positional QC, locus QC,
        and tNGS-specific QC rules as appropriate.

        Args:
            variants: List of Variant objects to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects
            target_coverage_map: Mapping of gene_id to TargetCoverage objects
            genes_with_valid_deletions: Mapping of gene_id to list of Variant objects with valid deletions

        Returns:
            The list of Variant objects with QC warnings applied
        """
        for variant in variants:
            logger.debug(f"Applying QC to {variant}")
            # Rule 4.2.1: Positional QC
            pos_qc_result = self._check_positional_qc(variant)
            self._update_variant_qc(variant, pos_qc_result)

            # Rule 4.2.2: Locus QC
            locus_qc_result = self._check_locus_qc(variant, pos_qc_result, locus_coverage_map, genes_with_valid_deletions)
            self._update_variant_qc(variant, locus_qc_result)

            # tNGS-specific QC (separate from rule structure)
            if self.config.TNGS:
                tngs_qc_result = self._check_tngs_qc(variant, locus_qc_result, locus_coverage_map)
                self._update_variant_qc(variant, tngs_qc_result)
        return variants

    def apply_wildtype_qc(self, variants: list[Variant], locus_coverage_map: Dict[str, LocusCoverage]) -> list[Variant]:
        """Rule 4.1: Applies QC to variants with no nucleotide change.

        Handles synthetic variants that represent unreported gene-drug associations.
        Also includes rule 4.2.2.3.1 (WT with insufficient coverage).

        Args:
            variants: List of Variant objects to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects

        Returns:
            The list of Variant objects with WT/NA QC applied
        """
        logger.debug(f"Applying WT QC results to {len(variants)} variants")
        for variant in variants:
            # initialize QC result with default None values
            qc_result = QCResult()

            locus_coverage = locus_coverage_map.get(variant.gene_id, None)

            # Handle NA/synthetic unreported variants
            if variant.nucleotide_change == "NA":
                qc_result.looker_interpretation="NA"
                qc_result.mdl_interpretation="NA"

                # WT Variant - if no nucleotide change, but exists in locus coverage map
                if locus_coverage:
                    qc_result.type = "WT"
                    qc_result.nucleotide_change = "WT"
                    qc_result.protein_change = "WT"

                    # if locus coverage is above threshold, set to WT interpretations
                    if not locus_coverage.has_breadth_below(self.config.MIN_PERCENT_COVERAGE):
                        qc_result.looker_interpretation = "S"
                        qc_result.mdl_interpretation = "WT"

                    # Rule 4.2.2.3.1: if locus coverage is below threshold, set to insufficient coverage
                    else:
                        qc_result.type = "Insufficient Coverage"
                        qc_result.looker_interpretation = "Insufficient Coverage"
                        qc_result.mdl_interpretation = "Insufficient Coverage"
                        qc_result.warning.add(self.LOCUS_QC_WARNING)
                        qc_result.fails_qc = True

            self._update_variant_qc(variant, qc_result)
        return variants

    def get_genes_with_valid_deletions(self, variants: list[Variant]) -> dict[str, list[Variant]]:
        """Helper function to get list of gene_ids with valid deletions.

        Args:
            variants: List of Variant objects to check

        Returns:
            Dictionary mapping gene_ids to lists of Variant objects with valid deletions
        """
        genes_with_valid_deletions = {}
        for variant in variants:
            if variant._is_valid_deletion():
                if variant.gene_id not in genes_with_valid_deletions:
                    genes_with_valid_deletions[variant.gene_id] = []
                genes_with_valid_deletions[variant.gene_id].append(variant)
        return genes_with_valid_deletions

    def _update_variant_qc(self, variant: Variant, qc_result: QCResult) -> None:
        """Apply a QCResult to a Variant object.

        Always sets:
            - `fails_qc`
            - `warning` (adds any warnings from the QCResult)

        Optionally, if set (not None), will overwrite:
            - `type`
            - `nucleotide_change`
            - `protein_change`
            - `rationale`
            - `confidence`
            - `looker_interpretation`
            - `mdl_interpretation`

        Args:
            variant: The variant to update
            qc_result: The QCResult to apply
        """
        # update/append any warnings and set/overwrite attributes if indicated by the QCResult
        qc_result.warning.update(variant.warning)
        updated_variant = variant.model_copy(
            update=qc_result.model_dump(exclude_none=True)
        )
        # Only log if there are changes to the variant (QC fail or warning added, or any attributes overwritten)
        changes = [
            attr for attr in qc_result.model_dump(exclude_none=True).keys()
            if (
                getattr(variant, attr) != getattr(updated_variant, attr) and
                getattr(variant, attr) is not None
            )
        ]
        qc_result_str = qc_result.repr_filtered([attr for attr in changes])
        if changes:
            logger.debug(f"Updating {variant} with QCResult: {qc_result_str}")

        for key, value in qc_result.model_dump(exclude_none=True).items():
            setattr(variant, key, value)

    def _check_positional_qc(self, variant: Variant) -> QCResult:
        """Rule 4.2.1 dispatcher: Route to deletion vs non-deletion positional QC.

        Args:
            variant: The variant to check

        Returns:
            QCResult with positional QC outcome
        """
        qc_result = QCResult(fails_qc=False)

        if not variant._is_valid_deletion():
            # rule 4.2.1.1 - Non-deletion positional QC
            if (
                variant.depth < self.config.MIN_DEPTH or
                variant.freq < self.config.MIN_FREQUENCY or
                variant.read_support < self.config.MIN_READ_SUPPORT
            ):
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS positional QC (rule 4.2.1.1)")

                qc_result.fails_qc = True
                qc_result.warning.add(self.POSITIONAL_QC_WARNING)
                return qc_result

            else:
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES positional QC")
                return qc_result

        # rule 4.2.1.2 - Deletion with some depth but below threshold
        elif variant.depth > 0 and variant.depth < self.config.MIN_DEPTH:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS positional QC (rule 4.2.1.2): Non-zero depth, but depth < {self.config.MIN_DEPTH}")
            qc_result.fails_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)
            return qc_result

        # rule 4.2.1.3 - Deletion with zero depth but passing frequency - TB Profiler quirk?
        elif variant.depth == 0 and variant.freq >= self.config.MIN_FREQUENCY:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES positional QC (rule 4.2.1.3): Zero depth, but freq >= {self.config.MIN_FREQUENCY}")
            return qc_result

        # Deletion passes positional QC (sufficient depth and frequency)
        logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES positional QC")
        return qc_result

    def _check_locus_qc(
        self,
        variant: Variant,
        qc_result: QCResult,
        locus_coverage_map: Dict[str, LocusCoverage],
        genes_with_valid_deletions: dict[str, list[Variant]],
    ) -> QCResult:
        """Rule 4.2.2 dispatcher: Check locus-level coverage QC.

        Args:
            variant: The variant to check
              qc_result: The result of the positional QC to inform locus QC decisions
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects
            genes_with_valid_deletions: List of gene_ids with valid deletions

        Returns:
            QCResult with locus QC outcome
        """

        locus_coverage = locus_coverage_map.get(variant.gene_id, None)
        if not locus_coverage:
            logger.debug(f"No locus coverage found for {variant.gene_name}|{variant.gene_id} at position {variant.pos}. Locus QC rules cannot be applied.")
            return qc_result

        # Check if gene contains variants with valid deletions
        gene_contains_valid_deletions = variant.gene_id in genes_with_valid_deletions
        if gene_contains_valid_deletions:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} contains valid deletion(s): {genes_with_valid_deletions[variant.gene_id]}")

        if locus_coverage.has_breadth_below(self.config.MIN_PERCENT_COVERAGE):
            # Rule 4.2.2.2: Low breadth of coverage with deletion present - PASS
            if gene_contains_valid_deletions:

                # NOTE: this conditional is in previous versions but NOT in the interpretation documentation. Not sure if intentional
                #  If this deletion also previously failed positional QC, add insufficient coverage warning?
                if variant.fails_qc:
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} contains valid deletion(s) with low breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}] and FAILS positional QC; FAILS locus QC; Adding `Insufficient Coverage` warning")
                    qc_result.fails_qc = True
                    qc_result.warning.add(self.LOCUS_QC_WARNING)
                    return qc_result

                # Rule 4.2.2.2: Low breadth of coverage with deletion present - PASS
                else:
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} contains valid deletion(s) with low breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}] and PASSES positional QC; PASSES locus QC (rule 4.2.2.2)")
                    return qc_result

            # Rule 4.2.2.3: Low breadth of coverage, no deletion
            # Note that Rule 4.2.2.3.1: WT/NA with low breadth of coverage is handled by `apply_wildtype_qc` function
            else:
                # Rule 4.2.2.3.2: S/U mutations with low breadth of coverage - FAIL
                if (
                    (variant.mdl_interpretation == "S") or
                    (variant.mdl_interpretation == "U")
                ):
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}]; FAILS locus QC; Adding `Insufficient Coverage` warning; Overwriting interpretation to `Insufficient Coverage` (rule 4.2.2.3.2)")
                    qc_result.fails_qc = True
                    qc_result.looker_interpretation = "Insufficient Coverage"
                    qc_result.mdl_interpretation = "Insufficient Coverage"
                    qc_result.warning.add(self.LOCUS_QC_WARNING)
                    return qc_result

                elif variant.mdl_interpretation == "R":
                    # Rule 4.2.2.3.3: R mutation with locus qc fail but NOT positional qc fail; add warning DO NOT not overwrite interpretation
                    if not variant.fails_qc:
                        logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}] and PASSES positional QC; Adding `Insufficient Coverage` warning (rule 4.2.2.3.3)")
                        qc_result.warning.add(self.LOCUS_QC_WARNING)
                        return qc_result
                    # Rule 4.2.2.3.4: R mutation with BOTH locus qc fail AND positional qc fail; add warning AND overwrite interpretation
                    else:
                        logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}] and FAILS positional QC; Adding `Insufficient Coverage` warning; Overwriting interpretation to `Insufficient Coverage` (rule 4.2.2.3.4)")
                        qc_result.fails_qc = True
                        qc_result.looker_interpretation = "Insufficient Coverage"
                        qc_result.mdl_interpretation = "Insufficient Coverage"
                        qc_result.warning.add(self.LOCUS_QC_WARNING)
                        return qc_result

                # not sure if this can happen
                else:
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}] but does not have a R/U/S mutation. Something went wrong. PASSES locus QC")
                    return qc_result
        else:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} has sufficient breadth of coverage [BC:{(locus_coverage.breadth_of_coverage):.3f}]; PASSES locus QC")
            return qc_result

    def _check_tngs_qc(
        self,
        variant: Variant,
        qc_result: QCResult,
        locus_coverage_map: Dict[str, LocusCoverage]
    ) -> QCResult:
        """Apply all tNGS-specific QC checks to a variant.

        Consolidates tNGS-specific QC, boundary checks, and outside-region checks.

        Args:
            variant: The variant to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects
        Returns:
            QCResult with tNGS QC outcome
        """
        # Check if mutation is outside tNGS primer regions
        locus_coverage = next(
            (lc for lc in locus_coverage_map.values()
            if variant.gene_id == lc.locus_tag and lc.contains_position(variant.pos)),
            None
        )
        if not locus_coverage:
            logger.debug(f"No tNGS locus coverage found for {variant.gene_name}|{variant.gene_id} at position {variant.pos}; FAILS tNGS QC; Adding warning: `This mutation is outside the expected region`")
            qc_result.fails_qc = True
            qc_result.warning.add("This mutation is outside the expected region")
            qc_result.rationale = "NA"
            qc_result.confidence = "NA"
            qc_result.looker_interpretation = "NA"
            qc_result.mdl_interpretation = "NA"

        # tNGS-specific gene/position checks
        if self._fails_tngs_specific_qc(variant):
            qc_result.fails_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)

        # tNGS boundary QC (read support vs frequency boundaries)
        if self._fails_tngs_boundary_qc(variant):
            qc_result.fails_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)

        return qc_result

    def _fails_tngs_specific_qc(self, variant: Variant) -> bool:
        """Check if a mutation (tNGS only) fails the tNGS-specific QC checks.

        Applies gene-specific frequency and read support thresholds for
        rrs, rrl, ethA (position 237), and rpoB (position 449).

        Args:
            variant: The variant to check

        Returns:
            True if the mutation fails QC, False if it passes
        """
        RRS_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRS_FREQUENCY"]
        RRS_READ_SUPPORT = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRS_READ_SUPPORT"]
        RRL_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRL_FREQUENCY"]
        RRL_READ_SUPPORT = self.config.TNGS_SPECIFIC_QC_OPTIONS["RRL_READ_SUPPORT"]
        ETHA237_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["ETHA237_FREQUENCY"]
        RPOB449_FREQUENCY = self.config.TNGS_SPECIFIC_QC_OPTIONS["RPOB449_FREQUENCY"]

        if variant.gene_name == "rrs":
            if variant.freq < RRS_FREQUENCY or variant.read_support < RRS_READ_SUPPORT:
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS tNGS-specific QC (rrs): freq < {RRS_FREQUENCY} or RS < {RRS_READ_SUPPORT}")
                return True
        elif variant.gene_name == "rrl":
            if variant.freq < RRL_FREQUENCY or variant.read_support < RRL_READ_SUPPORT:
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS tNGS-specific QC (rrl): freq < {RRL_FREQUENCY} or RS < {RRL_READ_SUPPORT}")
                return True
        elif variant.gene_name == "ethA" and Helper.get_position(variant.protein_change) == 237:
            if variant.freq < ETHA237_FREQUENCY:
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS tNGS-specific QC (ethA237): freq < {ETHA237_FREQUENCY}")
                return True
        elif variant.gene_name == "rpoB" and Helper.get_position(variant.protein_change) == 449:
            if variant.freq < RPOB449_FREQUENCY:
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS tNGS-specific QC (rpoB449): freq < {RPOB449_FREQUENCY}")
                return True
        logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES tNGS-specific QC")
        return False

    def _fails_tngs_boundary_qc(self, variant: Variant) -> bool:
        """Check if a mutation (tNGS only) fails the boundary QC checks.

        Applies read support and frequency boundary thresholds to determine
        if a variant passes tNGS QC.

        Args:
            variant: The variant to check

        Returns:
            True if the mutation fails QC, False if it passes
        """
        lower_rs = self.config.TNGS_READ_SUPPORT_BOUNDARIES[0]
        upper_rs = self.config.TNGS_READ_SUPPORT_BOUNDARIES[1]
        lower_f = self.config.TNGS_FREQUENCY_BOUNDARIES[0]
        upper_f = self.config.TNGS_FREQUENCY_BOUNDARIES[1]

        FAILS_QC = False
        if variant.read_support is None:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} is missing read support value; FAILS tNGS boundary QC: [RS:({lower_rs},{upper_rs}), F:({lower_f},{upper_f})]")
            return False

        if (lower_rs <= variant.read_support and variant.read_support < upper_rs):
            if (variant.freq < upper_f):
                FAILS_QC = True
        elif (variant.read_support >= upper_rs):
            if (variant.freq < lower_f):
                FAILS_QC = True
        elif (variant.read_support < lower_rs):
            FAILS_QC = True

        if FAILS_QC:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS tNGS boundary QC: [RS:({lower_rs},{upper_rs}), F:({lower_f},{upper_f})]")
        else:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES tNGS boundary QC: [RS:({lower_rs},{upper_rs}), F:({lower_f},{upper_f})]")
        return FAILS_QC
