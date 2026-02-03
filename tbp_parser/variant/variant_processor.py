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
        variants: List[Variant],
        sample_id: str,
        wildtype_candidates: List[str],
    ) -> List[Variant]:
        """Generate placeholder Variants for genes with no reported variants.

        Used for wild-type reporting - creates Variants for genes that
        had no mutations detected.

        Args:
            variants: List of Variant objects
        Returns:
            List of Variants for unreported genes
        """
        # function to determine WT/NA status. wildtype_candidates is a list of locus tags
        # found in the coverage map based on the user input BED file.
        wildtype_status = lambda gene_id: "WT" if gene_id in wildtype_candidates else "NA"

        all_unreported_variants = []
        unique_variants = set([variant.gene_id for variant in variants])
        unreported_variants = [v for v in GeneDatabase.GENE_DATABASE.keys() if v not in unique_variants]
        for var in unreported_variants:
            drug_list = GeneDatabase.get_drugs(var)
            for drug in drug_list:
                # setting attributes for a WT/NA (unreported) Variant
                gene_id = GeneDatabase.get_locus_tag(var)
                gene_name = GeneDatabase.get_gene_name(var)

                variant = Variant(
                    sample_id=sample_id,
                    pos=-1, # see `report_fmt` method in Variant class for handling of values during report writing
                    depth=-1,
                    freq=-1.0,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    type=wildtype_status(gene_id),
                    nucleotide_change=wildtype_status(gene_id),
                    protein_change=wildtype_status(gene_id),
                    confidence="NA",
                    drug=drug,
                    source="",
                    comment="",
                )
                # there's probably a better way of doing this. trying to keep the typing always consistent.
                # since read_support is freq * depth, to get it to report as NA in the lab report,
                # we need it to be negative. see `get_report_fmt` in lab_report.py
                variant.read_support = -1.0

                all_unreported_variants.append(variant)
        logger.debug(f"Generated {len(all_unreported_variants)} unreported variant-drug combinations for WT reporting")

        return all_unreported_variants

