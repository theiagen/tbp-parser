from dataclasses import dataclass
from typing import List

from tbp_parser.Variant.variant import Variant
from tbp_parser.Utilities.helper import Helper
import logging

logger = logging.getLogger(__name__)

@dataclass
class InterpretationResult:
    """Result of variant interpretation

    Attributes:
        confidence: WHO confidence classification or "No WHO annotation"
        looker_interpretation: Interpretation for Looker report (R, R-Interim, U, S, S-Interim)
        mdl_interpretation: Interpretation for MDL report (R, U, S, WT)
        rule_id: Document rule number for traceability (e.g., "1.1", "2.2.1.1")
        rationale: Human-readable description of the rule applied
    """
    confidence: str
    looker_interpretation: str
    mdl_interpretation: str
    rule_id: str
    rationale: str


class VariantInterpreter:
    """Class to handle interpretation of variants and expert rules"""

    RULE_TO_RATIONALE = {
        "1.2": "Expert rule 1.2. Novel drug targets",
        "2.2.1.1": "Expert rule 2.2.1.1. Loss-of-function",
        "2.2.2.1": "Expert rule 2.2.2.1. rpoB RRDR",
        "2.2.2.2": "Expert rule 2.2.2.2. rpoB non-RRDR",
        "3.2.1": "Expert rule 3.2.1. rrs",
        "3.2.2": "Expert rule 3.2.2. gyrA QRDR",
        "3.2.3": "Expert rule 3.2.3. gyrB QRDR",
        "3.2.4": "No WHO annotation or expert rule",
        "whov2": "Mutation in proximal promoter region",
    }

    def determine_interpretation(self, variants: List[Variant]) -> List[Variant]:
        """Updated interpretation method using rule-based structure.
        We only care about determining the interpretation for Variants of interest/sequenced with a nucleotide change.

        Args:
            variants: List of Variant objects to interpret
        Returns:
            List of Variant objects with updated interpretation fields
        """
        for variant in variants:
            # Handle synthetic unreported variants
            if not variant.nucleotide_change:
                variant.confidence = "NA"
                variant.looker_interpretation = "NA"
                variant.mdl_interpretation = "NA"
                variant.rationale = "NA"
                logger.debug(
                    f"{variant} has no mutation; marked as unreported variant; "
                    f"setting ['confidence', 'looker_interpretation', 'mdl_interpretation', 'rationale'] to 'NA'"
                )
                continue

            # Use interpret() method for all other cases
            result = self.interpret(variant)

            # Update Variant object with InterpretationResult attributes
            variant.confidence = result.confidence
            variant.looker_interpretation = result.looker_interpretation
            variant.mdl_interpretation = result.mdl_interpretation
            variant.rationale = result.rationale

            logger.debug(f"Final interpretation for {variant.gene_name} is: {result}")
        return variants

    def interpret(self, variant: Variant) -> InterpretationResult:
        """Main entry point for variant interpretation.

        Routes to appropriate section based on gene, mirroring interpretation document sections 1-3.

        Args:
            variant: The variant to interpret

        Returns:
            InterpretationResult with interpretation values and rule audit trail
        """
        # Section 1: CDC expert rule genes (new drugs)
        CDC_EXPERT_GENES = {"mmpR5", "Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC"}
        # Section 2: WHO expert rule genes
        WHO_EXPERT_GENES = {"katG", "pncA", "ethA", "gid", "rpoB"}

        gene = variant.gene_name
        logger.debug(f"Interpreting {variant}")

        # Section 1: CDC expert rule genes (mmpR5/Rv0678, atpE, pepQ, mmpL5, mmpS5, rrl, rplC)
        if gene in CDC_EXPERT_GENES:
            logger.debug(f"Gene {gene} matches CDC expert rules (Section 1)")
            return self._apply_section_1_rules(variant)

        # Section 2: WHO expert rule genes (katG, pncA, ethA, gid, rpoB)
        elif gene in WHO_EXPERT_GENES:
            logger.debug(f"Gene {gene} matches WHO expert rules (Section 2)")
            return self._apply_section_2_rules(variant)

        # Section 3: All other genes
        else:
            logger.debug(f"Gene {gene} uses other gene rules (Section 3)")
            return self._apply_section_3_rules(variant)

    # =========================================================================
    # SHARED HELPER METHODS (for rule consolidation)
    # =========================================================================

    def _has_who_confidence(self, variant: Variant) -> bool:
        """Check if variant has WHO confidence annotation.

        Returns True if confidence is not empty, not "Not found in WHO catalogue", and not "NA".
        """
        return (
            variant.confidence is not None and
            variant.confidence != "" and
            variant.confidence != "Not found in WHO catalogue" and
            variant.confidence != "No WHO annotation" and
            variant.confidence != "NA"
        )

    def _interpret_with_who_confidence(self, variant: Variant, rule_id: str) -> InterpretationResult:
        """Shared logic for rules 1.1, 2.1, 3.1 - when WHO confidence is available.

        Args:
            variant: The variant being interpreted
            rule_id: The rule ID for audit trail (e.g., "1.1", "2.1", "3.1")

        Returns:
            InterpretationResult with WHO confidence-based interpretation
        """
        # Maps WHO confidence values to (Looker interpretation, MDL interpretation)
        # Used by rules 1.1, 2.1, 3.1
        WHO_CONFIDENCE_MAPPING = {
            "Assoc w R": ("R", "R"),
            "Assoc w R - interim": ("R-Interim", "R"),
            "Assoc w R - Interim": ("R-Interim", "R"),
            "Uncertain significance": ("U", "U"),
            "Not assoc w R": ("S", "S"),
            "Not assoc w R - Interim": ("S-Interim", "S"),
            "Not assoc w R - interim": ("S-Interim", "S"),
        }

        confidence = variant.confidence
        if confidence in WHO_CONFIDENCE_MAPPING:
            looker, mdl = WHO_CONFIDENCE_MAPPING[confidence]
        else:
            # Fallback for unexpected confidence values
            logger.warning(f"Unknown WHO confidence value: {confidence}")
            looker, mdl = "X", "X"

        return InterpretationResult(
            confidence=confidence,
            looker_interpretation=looker,
            mdl_interpretation=mdl,
            rule_id=rule_id,
            rationale="WHO classification"
        )

    def _interpret_without_who_confidence(self, variant: Variant, rule_id: str) -> InterpretationResult:
        """Fallback for 'no WHO confidence' cases.

        Used by rules 1.2, 2.2.1.1 (non-LOF), 2.2.2.2, 3.2.4.

        Args:
            variant: The variant being interpreted
            rule_id: The rule ID for audit trail

        Returns:
            InterpretationResult based on mutation location and type
        """
        position_nt = Helper.get_position(variant.nucleotide_change)

        if variant._is_in_target_promoter(position_nt):
            logger.debug("Mutation in target promoter region")
            return InterpretationResult(
                confidence="No WHO annotation",
                looker_interpretation="U",
                mdl_interpretation="U",
                rule_id=rule_id,
                rationale=self.RULE_TO_RATIONALE.get('whov2', "MISSING RATIONALE") # proximal promoter rationale
            )
        elif variant._is_upstream_gene_variant():
            logger.debug("Mutation in upstream gene region")
            return InterpretationResult(
                confidence="No WHO annotation",
                looker_interpretation="S",
                mdl_interpretation="S",
                rule_id=rule_id,
                rationale=self.RULE_TO_RATIONALE.get(rule_id, "MISSING RATIONALE")
            )
        elif variant._is_in_orf() and variant._is_synonymous():
            logger.debug("Synonymous mutation in ORF")
            return InterpretationResult(
                confidence="No WHO annotation",
                looker_interpretation="S",
                mdl_interpretation="S",
                rule_id=rule_id,
                rationale=self.RULE_TO_RATIONALE.get(rule_id, "MISSING RATIONALE")
            )
        elif variant._is_in_orf() and not variant._is_synonymous():
            logger.debug("Nonsynonymous mutation in ORF")
            return InterpretationResult(
                confidence="No WHO annotation",
                looker_interpretation="U",
                mdl_interpretation="U",
                rule_id=rule_id,
                rationale=self.RULE_TO_RATIONALE.get(rule_id, "MISSING RATIONALE")
            )
        # not sure if this case is possible, but just in case
        else:
            return InterpretationResult(
                confidence="X",
                looker_interpretation="X",
                mdl_interpretation="X",
                rule_id=rule_id,
                rationale=self.RULE_TO_RATIONALE.get(rule_id, "MISSING RATIONALE")
            )

    # =========================================================================
    # SECTION 1 RULE METHODS - CDC Expert Rule Genes
    # =========================================================================

    def _apply_section_1_rules(self, variant: Variant) -> InterpretationResult:
        """Section 1: CDC genes (mmpR5/Rv0678, atpE, pepQ, mmpL5, mmpS5, rrl, rplC).

        These are genes related to new drugs based on CDC expert rules.
        """
        def _apply_rule_1_1_who_confidence(self, variant: Variant) -> InterpretationResult:
            """Rule 1.1: CDC genes WITH WHO confidence.

            Genes: mmpR5 (Rv0678), atpE, pepQ, mmpL5, mmpS5, rrl, rplC
            When WHO confidence is available, use that classification.
            """
            logger.debug(f"Applying rule 1.1 for {variant.gene_name} with WHO confidence")
            return self._interpret_with_who_confidence(variant, rule_id="1.1")

        def _apply_rule_1_2_no_who_confidence(self, variant: Variant) -> InterpretationResult:
            """Rule 1.2: CDC genes WITHOUT WHO confidence.

            Genes: mmpR5 (Rv0678), atpE, pepQ, mmpL5, mmpS5, rrl, rplC
            Interpretation based on location and mutation type.
            """

            logger.debug(f"Applying rule 1.2 for {variant.gene_name} without WHO confidence")

            def _apply_rule_1_2_rrl(self, variant: Variant) -> InterpretationResult:
                """Rule 1.2 special handling for rrl gene.

                rrl has specific nucleotide position ranges: 2003-2367 and 2449-3056
                """
                SPECIAL_POSITIONS = {
                    "rrl": [[2003, 2367], [2449, 3056]], # nt position ranges
                }

                position_nt = Helper.get_position(variant.nucleotide_change)

                # Check target promoter
                if variant._is_in_target_promoter(position_nt):
                    logger.debug("rrl mutation in target promoter; interpretation is 'U'")
                    return InterpretationResult(
                        confidence="No WHO annotation",
                        looker_interpretation="U",
                        mdl_interpretation="U",
                        rule_id="1.2",
                        rationale=self.RULE_TO_RATIONALE.get('whov2', "MISSING RATIONALE") # proximal promoter rationale
                    )

                # Check if in critical regions (nt positions 2003-2367 and 2449-3056)
                if Helper.is_mutation_within_range(position_nt, SPECIAL_POSITIONS["rrl"]):
                    logger.debug("rrl mutation in critical region (2003-2367 or 2449-3056); interpretation is 'U'")
                    return InterpretationResult(
                        confidence="No WHO annotation",
                        looker_interpretation="U",
                        mdl_interpretation="U",
                        rule_id="1.2",
                        rationale=self.RULE_TO_RATIONALE.get('1.2', "MISSING RATIONALE")
                    )

                # Outside critical regions and promoter
                logger.debug("rrl mutation outside critical regions and promoter; interpretation is 'S'")
                return InterpretationResult(
                    confidence="No WHO annotation",
                    looker_interpretation="S",
                    mdl_interpretation="S",
                    rule_id="1.2",
                    rationale=self.RULE_TO_RATIONALE.get('1.2', "MISSING RATIONALE")
                )

            # Special handling for rrl gene (has specific nucleotide position ranges)
            if variant.gene_name == "rrl":
                return _apply_rule_1_2_rrl(self, variant)

            # For other CDC genes: check promoter, upstream, then ORF
            logger.debug(f"Applying default interpretaion without WHO confidence for {variant.gene_name}")
            return self._interpret_without_who_confidence(variant, rule_id="1.2")

        if self._has_who_confidence(variant):
            return _apply_rule_1_1_who_confidence(self, variant)
        else:
            return _apply_rule_1_2_no_who_confidence(self, variant)

    # =========================================================================
    # SECTION 2 RULE METHODS - WHO Expert Rule Genes
    # =========================================================================

    def _apply_section_2_rules(self, variant: Variant) -> InterpretationResult:
        """Section 2: WHO genes (katG, pncA, ethA, gid, rpoB).

        These are genes covered by WHO expert rules.
        """
        def _apply_rule_2_1_who_confidence(self, variant: Variant) -> InterpretationResult:
            """Rule 2.1: WHO genes WITH WHO confidence.

            Genes: katG, pncA, ethA, gid, rpoB
            When WHO confidence is available, use that classification.
            """
            logger.debug(f"Applying rule 2.1 for {variant.gene_name} with WHO confidence")
            return self._interpret_with_who_confidence(variant, rule_id="2.1")

        def _apply_rule_2_2_no_who_confidence(self, variant: Variant) -> InterpretationResult:
            """Rule 2.2: WHO genes WITHOUT WHO confidence.

            Routes to loss-of-function (2.2.1) or rpoB (2.2.2) based on gene.
            """

            logger.debug(f"Applying rule 2.2 for {variant.gene_name} without WHO confidence")

            def _apply_rule_2_2_1_1_lof_in_orf(self, variant: Variant) -> InterpretationResult:
                """Rule 2.2.1.1: Loss-of-function expert rule for katG, pncA, ethA, gid, (NOT rpoB).

                If mutation is (loss-of-function) AND ((in ORF) or (within 30 bp upstream) -> R
                Otherwise, use default `interpret_without_who_confidence` interpretation.
                """
                logger.debug(f"Applying rule 2.2.1.1 (LOF check) for {variant.gene_name}")

                if (
                    (variant._is_loss_of_function()) and
                    (variant._is_in_orf() or variant._is_upstream_30_bp())
                ):
                    logger.debug("Loss-of-function mutation in ORF or within 30nt upstream; interpretation is 'R'")
                    return InterpretationResult(
                        confidence="No WHO annotation",
                        looker_interpretation="R",
                        mdl_interpretation="R",
                        rule_id="2.2.1.1",
                        rationale=self.RULE_TO_RATIONALE.get('2.2.1.1', "MISSING RATIONALE")
                    )

                # Fall through to default no-WHO interpretation - checks promoter, upstream, then ORF
                logger.debug(f"Not a qualifying LOF mutation for {variant.gene_name}")
                logger.debug(f"Applying default interpretaion without WHO confidence for {variant.gene_name}")
                return self._interpret_without_who_confidence(variant, rule_id="2.2.1.1")

            def _apply_rule_2_2_2_rpob(self, variant: Variant) -> InterpretationResult:
                """Rule 2.2.2: rpoB WITHOUT WHO confidence.

                Routes to RRDR (2.2.2.1) or non-RRDR (2.2.2.2) based on codon position.
                """

                logger.debug(f"Applying rule 2.2.2 for rpoB")

                def _apply_rule_2_2_2_1_rpob_rrdr(self, variant: Variant) -> InterpretationResult:
                    """Rule 2.2.2.1: rpoB mutation within RRDR (codons 426-452).

                    Nonsynonymous -> R, Synonymous -> S
                    """
                    logger.debug("Applying rule 2.2.2.1 for rpoB in RRDR region")

                    if variant._is_synonymous():
                        logger.debug("Synonymous mutation in RRDR; interpretation is 'S'")
                        return InterpretationResult(
                            confidence="No WHO annotation",
                            looker_interpretation="S",
                            mdl_interpretation="S",
                            rule_id="2.2.2.1",
                            rationale=self.RULE_TO_RATIONALE.get('2.2.2.1', "MISSING RATIONALE")
                        )
                    else:
                        logger.debug("Nonsynonymous mutation in RRDR; interpretation is 'R'")
                        return InterpretationResult(
                            confidence="No WHO annotation",
                            looker_interpretation="R",
                            mdl_interpretation="R",
                            rule_id="2.2.2.1",
                            rationale=self.RULE_TO_RATIONALE.get('2.2.2.1', "MISSING RATIONALE")
                        )

                def _apply_rule_2_2_2_2_rpob_non_rrdr(self, variant: Variant) -> InterpretationResult:
                    """Rule 2.2.2.2: rpoB mutation outside RRDR (codons 426-452).

                    Uses default no-WHO interpretation, but upstream variants → U (not S).
                    """
                    logger.debug("Applying rule 2.2.2.2 for rpoB outside RRDR; special interpretaion without WHO confidence")
                    return self._interpret_without_who_confidence(variant, rule_id="2.2.2.2")

                SPECIAL_POSITIONS = {
                    "rpoB": [426, 452], # RRDR codon range
                }
                position_aa = Helper.get_position(variant.protein_change)

                # Check if within RRDR (codons 426-452)
                if Helper.is_mutation_within_range(position_aa, SPECIAL_POSITIONS["rpoB"]):
                    return _apply_rule_2_2_2_1_rpob_rrdr(self, variant)
                else:
                    return _apply_rule_2_2_2_2_rpob_non_rrdr(self,variant)

            LOF_EXPERT_GENES = {"katG", "pncA", "ethA", "gid"}

            gene = variant.gene_name

            if gene in LOF_EXPERT_GENES:
                # note: skipping rule 2.2.1 because it only has one sub-rule (2.2.1.1)
                return _apply_rule_2_2_1_1_lof_in_orf(self, variant)
            elif gene == "rpoB":
                return _apply_rule_2_2_2_rpob(self, variant)
            else:
                logger.error(f"Gene {gene} not recognized in WHO expert rules 2.2")
                raise ValueError(f"Gene {gene} not recognized in WHO expert rules 2.2")

        if self._has_who_confidence(variant):
            return _apply_rule_2_1_who_confidence(self, variant)
        else:
            return _apply_rule_2_2_no_who_confidence(self, variant)

    # =========================================================================
    # SECTION 3 RULE METHODS - Other Genes
    # =========================================================================

    def _apply_section_3_rules(self, variant: Variant) -> InterpretationResult:
        """Section 3: All other genes.

        Includes rrs with special positions, gyrA/gyrB QRDR, and default rules.
        """
        def _apply_rule_3_1_who_confidence(self, variant: Variant) -> InterpretationResult:
            """Rule 3.1: Other genes WITH WHO confidence.

            When WHO confidence is available, use that classification.
            """
            logger.debug(f"Applying rule 3.1 for {variant.gene_name} with WHO confidence")
            return self._interpret_with_who_confidence(variant, rule_id="3.1")

        def _apply_rule_3_2_no_who_confidence(self, variant: Variant) -> InterpretationResult:
            """Rule 3.2: Other genes WITHOUT WHO confidence.

            Routes to specific gene rules or default fallback.
            """

            logger.debug(f"Applying rule 3.2 for {variant.gene_name}")

            def _apply_rule_3_2_1_rrs(self, variant: Variant) -> InterpretationResult:
                """Rule 3.2.1: rrs gene rules.

                Specific nucleotide positions 1401, 1402, 1484 → U
                Other positions → S
                """
                SPECIAL_POSITIONS = {
                    "rrs": [1401, 1402, 1484] # nt positions
                }

                logger.debug("Applying rule 3.2.1 for rrs")
                position_nt = Helper.get_position(variant.nucleotide_change)

                # Check critical positions (1401, 1402, 1484)
                critical_positions = SPECIAL_POSITIONS["rrs"]
                if any(pos in position_nt for pos in critical_positions):
                    logger.debug(f"rrs mutation at critical position; interpretation is 'U'")
                    return InterpretationResult(
                        confidence="No WHO annotation",
                        looker_interpretation="U",
                        mdl_interpretation="U",
                        rule_id="3.2.1",
                        rationale=self.RULE_TO_RATIONALE.get('3.2.1', "MISSING RATIONALE")
                    )
                # Other rrs positions
                else:
                    logger.debug("rrs mutation at non-critical position; interpretation is 'S'")
                    return InterpretationResult(
                        confidence="No WHO annotation",
                        looker_interpretation="S",
                        mdl_interpretation="S",
                        rule_id="3.2.1",
                        rationale=self.RULE_TO_RATIONALE.get('3.2.1', "MISSING RATIONALE")
                    )

            def _apply_rule_3_2_2_gyrA_qrdr(self, variant: Variant) -> InterpretationResult:
                """Rule 3.2.2: gyrA QRDR rules.

                Codons 88-94, nonsynonymous → U
                Otherwise use default.
                """
                SPECIAL_POSITIONS = {
                    "gyrA": [88, 94], # QRDR codon range
                }
                logger.debug("Applying rule 3.2.2 for gyrA")
                position_aa = Helper.get_position(variant.protein_change)

                # Check if in QRDR (codons 88-94) and nonsynonymous
                if Helper.is_mutation_within_range(position_aa, SPECIAL_POSITIONS["gyrA"]):
                    if not variant._is_synonymous():
                        logger.debug("gyrA nonsynonymous mutation in QRDR; interpretation is 'U'")
                        return InterpretationResult(
                            confidence="No WHO annotation",
                            looker_interpretation="U",
                            mdl_interpretation="U",
                            rule_id="3.2.2",
                            rationale=self.RULE_TO_RATIONALE.get('3.2.2', "MISSING RATIONALE")
                        )

                # Fall through to default
                return _apply_rule_3_2_4_no_who_confidence(self, variant)

            def _apply_rule_3_2_3_gyrB_qrdr(self, variant: Variant) -> InterpretationResult:
                """Rule 3.2.3: gyrB QRDR rules.

                Codons 446-507, nonsynonymous → U
                Otherwise use default.
                """
                SPECIAL_POSITIONS = {
                    "gyrB": [446, 507],
                }

                logger.debug("Applying rule 3.2.3 for gyrB")
                position_aa = Helper.get_position(variant.protein_change)

                # Check if in QRDR (codons 446-507) and nonsynonymous
                if Helper.is_mutation_within_range(position_aa, SPECIAL_POSITIONS["gyrB"]):
                    if not variant._is_synonymous():
                        logger.debug("gyrB nonsynonymous mutation in QRDR; interpretation is 'U'")
                        return InterpretationResult(
                            confidence="No WHO annotation",
                            looker_interpretation="U",
                            mdl_interpretation="U",
                            rule_id="3.2.3",
                            rationale=self.RULE_TO_RATIONALE.get('3.2.3', "MISSING RATIONALE")
                        )

                # Fall through to default
                return _apply_rule_3_2_4_no_who_confidence(self, variant)

            def _apply_rule_3_2_4_no_who_confidence(self, variant: Variant) -> InterpretationResult:
                """Rule 3.2.4: Default fallback for other genes without WHO confidence.

                Uses standard no-WHO interpretation based on mutation type and location.
                """
                logger.debug(f"Applying default interpretaion without WHO confidence for {variant.gene_name}")
                return self._interpret_without_who_confidence(variant, rule_id="3.2.4")

            gene = variant.gene_name

            # rrs has specific position rules
            if gene == "rrs" and variant._is_in_orf():
                return _apply_rule_3_2_1_rrs(self, variant)

            # gyrA has QRDR rules
            if gene == "gyrA" and variant._is_in_orf():
                return _apply_rule_3_2_2_gyrA_qrdr(self, variant)
            # gyrB has QRDR rules
            if gene == "gyrB" and variant._is_in_orf():
                return _apply_rule_3_2_3_gyrB_qrdr(self, variant)

            # Default fallback for all other genes `interpret_without_who_confidence`
            return _apply_rule_3_2_4_no_who_confidence(self, variant)

        if self._has_who_confidence(variant):
            return _apply_rule_3_1_who_confidence(self, variant)
        else:
            return _apply_rule_3_2_no_who_confidence(self, variant)
