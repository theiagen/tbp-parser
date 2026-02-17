
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

    def apply_qc(self, variants: list[Variant], locus_coverage_map: Dict[str, LocusCoverage]) -> list[Variant]:
        """Apply all QC checks to a list of variants.

        Iterates through variants and applies positional QC, locus QC,
        and tNGS-specific QC rules as appropriate.

        Args:
            variants: List of Variant objects to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects

        Returns:
            The list of Variant objects with QC warnings applied
        """
        for variant in variants:
            # Rule 4.2.1: Positional QC
            pos_qc_result = self._check_positional_qc(variant)
            self._update_variant_qc(variant, pos_qc_result)

            # Rule 4.2.2: Locus QC
            locus_qc_result = self._check_locus_qc(variant, pos_qc_result, locus_coverage_map)
            self._update_variant_qc(variant, locus_qc_result)

            # tNGS-specific QC (separate from rule structure)
            if self.config.TNGS:
                self._apply_tngs_qc(variant, locus_coverage_map)

            # Track valid deletions (passed all QC, is a deletion)
            if not variant.fails_qc and self._is_deletion(variant):
                setattr(variant, "is_valid_deletion", True)

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
        # initialize QC result with default None values

        for variant in variants:
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

    def _update_variant_qc(self, variant: Variant, qc_result: QCResult) -> None:
        """Apply a QCResult to a Variant object.

        Always sets:
            - `fails_qc`
            - `warning` (adds any warnings from the QCResult)

        Optionally, if set (not None), will overwrite:
            - `type`
            - `nucleotide_change`
            - `protein_change`
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

    def _is_deletion(self, variant: Variant) -> bool:
        """Check if a variant represents a deletion.

        Args:
            variant: The variant to check

        Returns:
            True if the nucleotide_change contains 'del'
        """
        return "del" in variant.nucleotide_change

    def _check_positional_qc(self, variant: Variant) -> QCResult:
        """Rule 4.2.1 dispatcher: Route to deletion vs non-deletion positional QC.

        Args:
            variant: The variant to check

        Returns:
            QCResult with positional QC outcome
        """
        qc_result = QCResult(fails_qc=False)

        if not self._is_deletion(variant):
            # rule 4.2.1.1 - Non-deletion positional QC
            if (
                variant.depth < self.config.MIN_DEPTH or
                variant.freq < self.config.MIN_FREQUENCY or
                variant.read_support < self.config.MIN_READ_SUPPORT
            ):
                logger.debug(f"{variant} fails positional QC (rule 4.2.1.1): depth={variant.depth}, freq={variant.freq}, RS={variant.read_support}")
                qc_result.fails_qc = True
                qc_result.warning.add(self.POSITIONAL_QC_WARNING)
                return qc_result

        # rule 4.2.1.2 - Deletion with some depth but below threshold
        elif variant.depth > 0 and variant.depth < self.config.MIN_DEPTH:
            logger.debug(f"{variant} deletion fails positional QC (rule 4.2.1.2): depth={variant.depth} > 0 but < {self.config.MIN_DEPTH}")
            qc_result.fails_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)
            return qc_result

        # rule 4.2.1.3 - Deletion with zero depth but passing frequency - TB Profiler quirk
        elif variant.depth == 0 and variant.freq >= self.config.MIN_FREQUENCY:
            logger.debug(f"{variant} deletion passes positional QC (rule 4.2.1.3): depth=0, freq={variant.freq} >= {self.config.MIN_FREQUENCY}")
            return qc_result

        # Deletion passes positional QC (sufficient depth and frequency)
        return qc_result

    def _check_locus_qc(self, variant: Variant, qc_result: QCResult, locus_coverage_map: Dict[str, LocusCoverage]) -> QCResult:
        """Rule 4.2.2 dispatcher: Check locus-level coverage QC.

        Args:
            variant: The variant to check
              qc_result: The result of the positional QC to inform locus QC decisions
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects

        Returns:
            QCResult with locus QC outcome
        """
        locus_coverage = locus_coverage_map.get(variant.gene_id, None)

        if locus_coverage and locus_coverage.has_breadth_below(self.config.MIN_PERCENT_COVERAGE):
            # Rule 4.2.2.2: Low breadth of coverage with deletion present - PASS
            if self._is_deletion(variant):

                # NOTE: this conditional is in previous versions but NOT in the interpretation documentation. Not sure if intentional
                #  If this deletion also previously failed positional QC, add insufficient coverage warning?
                if variant.fails_qc:
                    logger.debug(f"{variant} deletion with low breadth AND positional fail. Adding insufficient coverage warning.")
                    qc_result.fails_qc = True
                    qc_result.warning.add(self.LOCUS_QC_WARNING)
                    return qc_result

                # Rule 4.2.2.2: Low breadth of coverage with deletion present - PASS
                else:
                    logger.debug(f"{variant} deletion with low breadth but no positional fail. Passes locus QC.")
                    return qc_result

            # Rule 4.2.2.3: Low breadth of coverage, no deletion
            # Note that Rule 4.2.2.3.1: WT/NA with low breadth of coverage is handled by `apply_wildtype_qc` function
            else:
                # Rule 4.2.2.3.2: S/U or R (QC_RESISTANT_MUTATIONS) with low breadth of coverage - FAIL
                if variant.looker_interpretation == "R" and self.config.QC_RESISTANT_MUTATIONS and not variant.fails_qc:
                    logger.debug(f"{variant} non-deletion with low breadth of coverage but R interpretation did NOT fail positional QC. Only adding insufficient coverage warning.")
                    qc_result.warning.add(self.LOCUS_QC_WARNING)
                    return qc_result

                elif (
                    (variant.mdl_interpretation == "S") or
                    (variant.mdl_interpretation == "U") or
                    (variant.mdl_interpretation == "R" and self.config.QC_RESISTANT_MUTATIONS and variant.fails_qc)
                ):
                    logger.debug(f"{variant} non-deletion with low breadth of coverage. Adding insufficient coverage warning and overwriting interpretations.")
                    qc_result.fails_qc = True
                    qc_result.looker_interpretation = "Insufficient Coverage"
                    qc_result.mdl_interpretation = "Insufficient Coverage"
                    qc_result.warning.add(self.LOCUS_QC_WARNING)
                    return qc_result

                else:
                    return qc_result
        else:
            return qc_result

    def _apply_tngs_qc(self, variant: Variant, locus_coverage_map: Dict[str, LocusCoverage]) -> None:
        """Apply all tNGS-specific QC checks to a variant.

        Consolidates tNGS-specific QC, boundary checks, and outside-region checks.
        Modifies variant directly (does not use QCResult pattern).

        Args:
            variant: The variant to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects
        """
        # tNGS-specific gene/position checks
        if self._fails_tngs_specific_qc(variant):
            variant.fails_qc = True
            variant.warning.add(self.POSITIONAL_QC_WARNING)

        # tNGS boundary QC (read support vs frequency boundaries)
        if self._fails_tngs_boundary_qc(variant):
            variant.fails_qc = True
            variant.warning.add(self.POSITIONAL_QC_WARNING)

        # Check if mutation is outside tNGS primer regions
        locus_coverage = locus_coverage_map.get(variant.gene_id, None)
        if (
            not locus_coverage or
            not locus_coverage.contains_position(variant.pos)
        ):
            variant.fails_qc = True
            variant.warning.add("This mutation is outside the expected region")
            variant.rationale = "NA"
            variant.confidence = "NA"
            variant.looker_interpretation = "NA"
            variant.mdl_interpretation = "NA"

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
