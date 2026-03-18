from __future__ import annotations
from pydantic import BaseModel, Field
from typing import Dict, Optional, TYPE_CHECKING
from tbp_parser.Utilities.config import Configuration
import logging

if TYPE_CHECKING:
    from tbp_parser.Coverage.coverage_data import LocusCoverage, TargetCoverage
    from tbp_parser.Variant.variant import Variant


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
    fails_positional_qc: bool = False
    fails_locus_qc: bool = False
    warning: set[str] = Field(default_factory=set)

    def repr_filtered(self, attr_list: list[str]) -> str:
        """Return a string representation of the QCResult, filtered to only include specified attributes."""
        attrs = {k: v for k, v in self.model_dump(exclude_none=True).items() if k in attr_list}
        return f"QCResult({', '.join([f'{k}={v}' for k, v in attrs.items()])})"

class VariantQC:

    POSITIONAL_QC_WARNING = "Failed quality in the mutation position"
    LOCUS_QC_WARNING = "Insufficient coverage in locus"

    def __init__(self):
        self.config = Configuration.get_instance()

    def qc(
        self,
        variants: list[Variant],
        unreported_variants: list[Variant],
        locus_coverage_map: Dict[str, LocusCoverage],
        target_coverage_map: Dict[str, TargetCoverage],
    ) -> list[Variant]:
        """
        Main method for performing QC on a list of Variant objects.
        Applies positional QC, locus QC, and tNGS-specific QC rules to all variants.
        Also applies WT/NA QC to unreported variants and identifies genes with valid deletions.
        Args:
            variants: List of Variant objects to QC
            unreported_variants: List of unreported Variant objects to QC (WT/NA QC only)
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects for locus QC
            target_coverage_map: Mapping of gene_name to TargetCoverage objects for assigning valid deletions
        Returns:
            List of all Variant objects with QC applied
        """
        # Apply QC to all variants. Note only LocusCoverage is used for QC
        variants = self.apply_qc(variants, locus_coverage_map, target_coverage_map)

        # Apply WT/NA QC to unreported_variants
        unreported_variants = self.apply_wildtype_qc(unreported_variants, locus_coverage_map)

        # merge variants and unreported_variants for reporting
        reported_variants = variants + unreported_variants
        return reported_variants

    def assign_variants_with_valid_deletions(
        self,
        variant: Variant,
        locus_coverage_map: Dict[str, LocusCoverage],
        target_coverage_map: Dict[str, TargetCoverage],
    ) -> None:
        """Helper function to assign list of Variants with valid deletions that fall within a specific coverage region.

        Args:
            variant: Variant object to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects
            target_coverage_map: Mapping of gene_name to TargetCoverage objects
        """
        # If not a deletion or the variant fails positional QC, skip it
        if not variant._is_deletion_in_orf() or variant.fails_positional_qc:
            return

        # Determine absolute start and end positions of the variant for coverage overlap checks.
        var_abs_start = variant.absolute_start if isinstance(variant.absolute_start, int) else variant.pos
        var_abs_end = variant.absolute_end if isinstance(variant.absolute_end, int) else variant.pos

        # Map each coverage_map to the variant attribute used as its key
        for coverage_map in [target_coverage_map, locus_coverage_map]:
            for coverage in coverage_map.values():
                # Check if the variant is a valid deletion and falls within the coverage region
                if (
                    variant.gene_id == coverage.locus_tag and
                    coverage.overlaps_range(var_abs_start, var_abs_end)
                ):
                    coverage.valid_deletions.append(variant)
                    logger.debug(f"Assigned {variant.gene_name}|{variant.gene_id} as a valid deletion within the coverage region of the {coverage.__class__.__name__} object")
                    # Check if ERR coverage exists and assign Variants to valid_deletions if position falls within ERR region
                    if coverage.err_coverage and coverage.err_coverage.overlaps_range(var_abs_start, var_abs_end):
                        coverage.err_coverage.valid_deletions.append(variant)
                        logger.debug(f"Assigned {variant.gene_name}|{variant.gene_id} as a valid deletion within the ERR coverage region of the {coverage.err_coverage.__class__.__name__} object")
        return

    def apply_qc(
        self,
        variants: list[Variant],
        locus_coverage_map: Dict[str, LocusCoverage],
        target_coverage_map: Dict[str, TargetCoverage],
    ) -> list[Variant]:
        """Apply all QC checks to a list of variants.

        Iterates through variants and applies positional QC, locus QC,
        and tNGS-specific QC rules as appropriate.

        Args:
            variants: List of Variant objects to check
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects
            target_coverage_map: Mapping of gene_name to TargetCoverage objects

        Returns:
            The list of Variant objects with QC warnings applied
        """
        for variant in variants:
            logger.debug(f"Applying QC to {variant}")

            # Rule 4.2.1: Positional QC
            qc_result = self.check_positional_qc(variant)
            variant = self.update_variant_qc_result(variant, qc_result)

            # tNGS-specific QC (separate from rule structure)
            if self.config.TNGS:
                tngs_qc_result = self.check_tngs_qc(variant, qc_result, locus_coverage_map)
                variant = self.update_variant_qc_result(variant, tngs_qc_result)

            # Find/assign valid deletions to coverage objects if the variant passes positional QC and is within the coverage region
            self.assign_variants_with_valid_deletions(variant, locus_coverage_map, target_coverage_map)

            # Rule 4.2.2: Locus QC
            qc_result = self.check_locus_qc(variant, qc_result, locus_coverage_map)
            variant = self.update_variant_qc_result(variant, qc_result)

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
                        qc_result.fails_locus_qc = True

            variant = self.update_variant_qc_result(variant, qc_result)
        return variants

    def update_variant_qc_result(self, variant: Variant, qc_result: QCResult) -> Variant:
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
        return variant

    def check_positional_qc(self, variant: Variant) -> QCResult:
        """Rule 4.2.1 dispatcher: Route to deletion vs non-deletion positional QC.

        Args:
            variant: The variant to check

        Returns:
            QCResult with positional QC outcome
        """
        qc_result = QCResult(fails_positional_qc=False)

        if not variant._is_deletion_in_orf():
            # rule 4.2.1.1 - Non-deletion positional QC
            if (
                variant.depth < self.config.MIN_DEPTH or
                variant.freq < self.config.MIN_FREQUENCY or
                variant.read_support < self.config.MIN_READ_SUPPORT
            ):
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS positional QC (rule 4.2.1.1)")

                qc_result.fails_positional_qc = True
                qc_result.warning.add(self.POSITIONAL_QC_WARNING)
                return qc_result

            else:
                logger.debug(f"{variant.gene_name}|{variant.gene_id} [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES positional QC")
                return qc_result

        # rule 4.2.1.2 - Deletion with some depth but below threshold
        elif variant.depth > 0 and variant.depth < self.config.MIN_DEPTH:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS positional QC (rule 4.2.1.2): Non-zero depth, but depth < {self.config.MIN_DEPTH}")
            qc_result.fails_positional_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)
            return qc_result

        # rule 4.2.1.3 - Deletion with zero depth but passing frequency - TB Profiler quirk?
        elif variant.depth == 0 and variant.freq >= self.config.MIN_FREQUENCY:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES positional QC (rule 4.2.1.3): Zero depth, but freq >= {self.config.MIN_FREQUENCY}")
            return qc_result

        # NOTE: not in interpretation document. But i think this shouldn't pass positional QC if depth is 0 and freq is below threshold
        elif variant.depth == 0 and variant.freq < self.config.MIN_FREQUENCY:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; FAILS positional QC (rule 4.2.1.4): Zero depth and freq < {self.config.MIN_FREQUENCY}")
            qc_result.fails_positional_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)
            return qc_result

        # Deletion passes positional QC (sufficient depth and frequency)
        logger.debug(f"{variant.gene_name}|{variant.gene_id} is a valid deletion with [D:{variant.depth}, RS:{(variant.read_support):.0f}, F:{(variant.freq):.3f}]; PASSES positional QC")
        return qc_result

    def check_locus_qc(
        self,
        variant: Variant,
        qc_result: QCResult,
        locus_coverage_map: Dict[str, LocusCoverage],
    ) -> QCResult:
        """Rule 4.2.2 dispatcher: Check locus-level coverage QC.

        Args:
            variant: The variant to check
            qc_result: The result of the positional QC to inform locus QC decisions
            locus_coverage_map: Mapping of gene_id to LocusCoverage objects

        Returns:
            QCResult with locus QC outcome
        """
        locus_coverage = locus_coverage_map.get(variant.gene_id, None)
        if not locus_coverage:
            logger.debug(f"No locus coverage found for {variant.gene_name}|{variant.gene_id} at position {variant.pos}. Locus QC rules cannot be applied.")
            return qc_result

        # if ERR coverage exists AND the `--use_err_as_brr` flag is set, use ERR coverage for locus QC instead of overall locus coverage
        if self.config.USE_ERR_AS_BRR and locus_coverage.err_coverage:
            logger.debug(f"Using ERR coverage for locus QC of {variant.gene_name}|{variant.gene_id}")
            locus_coverage = locus_coverage.err_coverage

        # Determine if locus coverage is below threshold for QC fail and check if variant is a valid deletion
        has_low_boc = locus_coverage.has_breadth_below(self.config.MIN_PERCENT_COVERAGE)
        boc = locus_coverage.breadth_of_coverage
        has_valid_deletion = locus_coverage.contains_loci_with_valid_deletion(variant.gene_id)

        if has_low_boc:
            # Always add warning if locus has low breadth of coverage
            qc_result.warning.add(self.LOCUS_QC_WARNING)

            # Rule 4.2.2.2: Low breadth of coverage with deletion present - PASS
            if has_valid_deletion:
                # NOTE: this conditional is in previous versions but NOT in the interpretation documentation. Not sure if intentional
                #  If this deletion also previously failed positional QC, add insufficient coverage warning?
                if variant.fails_positional_qc:
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} contains valid deletion(s) with low breadth of coverage [BC:{(boc):.3f}] and FAILS positional QC; FAILS locus QC; Adding `Insufficient Coverage` warning")
                    qc_result.fails_locus_qc = True
                    return qc_result

                # Rule 4.2.2.2: Low breadth of coverage with deletion present - PASS
                else:
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} contains valid deletion(s) with low breadth of coverage [BC:{(boc):.3f}] and PASSES positional QC; PASSES locus QC (rule 4.2.2.2)")
                    return qc_result

            # Rule 4.2.2.3: Low breadth of coverage, no deletion
            # Note that Rule 4.2.2.3.1: WT/NA with low breadth of coverage is handled by `apply_wildtype_qc` function
            else:
                # Rule 4.2.2.3.2: S/U mutations with low breadth of coverage - FAIL
                if (variant.mdl_interpretation == "S") or (variant.mdl_interpretation == "U"):
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(boc):.3f}]; FAILS locus QC; Adding `Insufficient Coverage` warning; Overwriting interpretation to `Insufficient Coverage` (rule 4.2.2.3.2)")
                    qc_result.fails_locus_qc = True
                    qc_result.looker_interpretation = "Insufficient Coverage"
                    qc_result.mdl_interpretation = "Insufficient Coverage"
                    return qc_result

                elif variant.mdl_interpretation == "R":
                    # Rule 4.2.2.3.3: R mutation with locus qc fail but NOT positional qc fail; add warning DO NOT not overwrite interpretation
                    if not variant.fails_positional_qc:
                        logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(boc):.3f}] and PASSES positional QC; Adding `Insufficient Coverage` warning (rule 4.2.2.3.3)")
                        return qc_result

                    # Rule 4.2.2.3.4: R mutation with BOTH locus qc fail AND positional qc fail; add warning AND overwrite interpretation
                    else:
                        logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(boc):.3f}] and FAILS positional QC; Adding `Insufficient Coverage` warning; Overwriting interpretation to `Insufficient Coverage` (rule 4.2.2.3.4)")
                        qc_result.fails_locus_qc = True
                        qc_result.looker_interpretation = "Insufficient Coverage"
                        qc_result.mdl_interpretation = "Insufficient Coverage"
                        return qc_result

                # not sure if this can happen
                else:
                    logger.debug(f"{variant.gene_name}|{variant.gene_id} has low breadth of coverage [BC:{(boc):.3f}] but does not have a R/U/S mutation. Something went wrong. PASSES locus QC")
                    return qc_result
        else:
            logger.debug(f"{variant.gene_name}|{variant.gene_id} has sufficient breadth of coverage [BC:{(boc):.3f}]; PASSES locus QC")
            return qc_result

    def check_tngs_qc(
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
            if variant.gene_id == lc.locus_tag and lc.overlaps_range(
                variant.absolute_start if isinstance(variant.absolute_start, int) else variant.pos,
                variant.absolute_end if isinstance(variant.absolute_end, int) else variant.pos,
            )),
            None
        )
        if not locus_coverage:
            logger.debug(f"No tNGS locus coverage found for {variant.gene_name}|{variant.gene_id} at position {variant.pos}; FAILS tNGS QC; Adding warning: `This mutation is outside the expected region`")
            qc_result.fails_locus_qc = True
            qc_result.warning.add("This mutation is outside the expected region")
            qc_result.rationale = "NA"
            qc_result.confidence = "NA"
            qc_result.looker_interpretation = "NA"
            qc_result.mdl_interpretation = "NA"

        # tNGS boundary QC (read support vs frequency boundaries)
        if self._fails_tngs_boundary_qc(variant):
            qc_result.fails_positional_qc = True
            qc_result.warning.add(self.POSITIONAL_QC_WARNING)

        return qc_result

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
