import logging
from typing import List

from variant import Variant
from utils.gene_database import GeneDatabase

logger = logging.getLogger(__name__)

class VariantProcessor:
    """Process and analyze variants."""

    @staticmethod
    def deduplicate_variants(variants: List[Variant]) -> List[Variant]:
        """Deduplicate variants, keeping the one with best annotation.

        When multiple identical variants exist (same gene, position, drug),
        keeps the one with the highest WHO confidence ranking.

        Args:
            variants: List of Variant objects (may contain duplicates)

        Returns:
            List of deduplicated Variants
        """
        unique_variants = {}

        logger.debug(f"Deduplicating Variant objects based on `gene_id`, `gene_name`, `type`, `nucleotide_change`, and `drug`")
        for variant in variants:
            # create unique key based on identifying attributes
            key = (
                variant.gene_id,
                variant.gene_name,
                variant.type,
                variant.nucleotide_change,
                variant.drug
            )

            if key not in unique_variants:
                unique_variants[key] = variant
            elif variant.is_better_annotation_than(unique_variants[key]):
                unique_variants[key] = variant

        result = list(unique_variants.values())
        logger.debug(
            f"Deduplicated {len(variants)} variants to {len(result)} unique variants"
        )

        return result

    @staticmethod
    def generate_unreported_variants(
        variants: List[Variant]
    ) -> List[Variant]:
        """Generate placeholder Variants for genes with no reported variants.

        Used for wild-type reporting - creates Variants for genes that
        had no mutations detected.

        Args:
            variants: List of Variant objects (may contain duplicates)
        Returns:
            List of placeholder Variants for unreported genes
        """
        all_unreported_variants = []
        unique_variants = set([variant.gene_id for variant in variants])
        unreported_variants = [v for v in GeneDatabase.GENE_DATABASE.keys() if v not in unique_variants]
        for var in unreported_variants:
            for drug in GeneDatabase.get_drugs(var):
                variant = Variant.from_thin_air(
                    gene_id=GeneDatabase.get_locus_tag(var),
                    gene_name=GeneDatabase.get_gene_name(var),
                    drug=drug,
                )
                all_unreported_variants.append(variant)
        logger.debug(f"Generated {len(all_unreported_variants)} unreported variant-drug combinations for WT reporting")

        return variants + all_unreported_variants

