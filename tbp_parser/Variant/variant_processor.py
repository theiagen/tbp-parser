import logging
from typing import List, Optional, Tuple

from Variant import Variant, VariantRecord, Annotation
from Utilities import GeneDatabase

logger = logging.getLogger(__name__)

class VariantProcessor:
    """Process VariantRecords into Variant objects."""

    def process(self, variant_records: List[VariantRecord], sample_id: str) -> Tuple[List[Variant], List[Variant]]:
        """Main method for processing VariantRecords into Variant objects.

        Args:
            variant_records: List of VariantRecord objects
            sample_id: Sample identifier for generating unreported variants
        Returns:
            Tuple of (processed Variants, unreported Variants)
        """
        # Step 1: Process VariantRecords into Variants (expansion + extraction)
        reported_variants = self.process_variant_records(variant_records)

        # Step 2: Deduplicate Variants, keeping the one with best annotation
        reported_variants = self.deduplicate_variants(reported_variants)

        # Step 3: Generate unreported_variants for gene/drug associations not in original set
        unreported_variants = self.generate_unreported_variants(reported_variants, sample_id)
        return reported_variants, unreported_variants

    def process_variant_records(self, variant_records: List[VariantRecord]) -> List[Variant]:
        """Process VariantRecords into Variant objects.

        Expands consequences for each variant entry,
        and extracts individual Variant objects (one per drug annotation).

        Args:
            variant_records: List of VariantRecord objects
        Returns:
            List of Variant objects extracted from the VariantRecords
        """
        reported_variants = []
        for variant_record in variant_records:
            # expand consequences in VariantRecords if applicable -> will always return list of VariantRecords
            expanded_variant_records = self._expand_consequences(variant_record)

            # extract Variant objects from each expanded VariantRecord based on annotations -> will always return list of Variant objects
            for record in expanded_variant_records:
                variants = self._get_variants_from_annotations(record)
                reported_variants.extend(variants)

        logger.debug(f"Total Variant objects from VariantRecords: {len(reported_variants)}")
        return reported_variants

    def _expand_consequences(self, variant_record: VariantRecord) -> List[VariantRecord]:
        """Expand a VariantRecord into multiple entries based on consequences.

        Some variants affect multiple genes (e.g., a mutation in mmpR5 that
        also affects mmpL5 and mmpS5). This creates separate entries for each
        affected gene. Always returns at least the original VariantRecord, even if no valid consequences exist.

        Args:
            variant_record: Single VariantRecord entry from TBProfiler JSON
        Returns:
            List of VariantRecords (original + consequences)
        """
        # Only expand if consequences exist and gene is one of mmpR5/mmpL5/mmpS5 otherwise return original
        all_variant_records = [variant_record]
        if (
            not variant_record.consequences or
            variant_record.gene_id not in ['Rv0676c', 'Rv0677c', 'Rv0678']
        ):
            return all_variant_records

        logger.debug(f"Expanding {len(variant_record.consequences)} consequences for VariantRecord entry: {str(variant_record)}")

        for consequence in variant_record.consequences:
            # skip and do not expand consequence if it's for the same gene as the original VariantRecord
            if consequence.gene_id == variant_record.gene_id:
                continue
            # Create a copy and update with consequence-specific data
            new_vr = VariantRecord.from_consequences(variant_record, consequence)
            all_variant_records.append(new_vr)

        return all_variant_records

    def _get_variants_from_annotations(self, variant_record: 'VariantRecord') -> List['Variant']:
        """Creates Variants from a VariantRecord, expanding to all relevant drug associations
        Args:
            variant_record: A dictionary representing a single variant record from the TBProfiler JSON output.
        Returns:
            A list of Variant instances.
        """
        if not variant_record.annotation:
            return []

        # Build complete set of annotations (original + synthetic)
        all_annotations = self._expand_annotations_for_all_drugs(
            variant_record.annotation,
            variant_record.gene_associated_drugs,
            variant_record.gene_id
        )
        # Convert each annotation to a Variant
        variants = [
            Variant(**variant_record.model_dump(), **anno.model_dump())
            for anno in all_annotations
        ]

        logger.debug(
            f"Created {len(variants)} Variant objects from {len(variant_record.annotation)} "
            f"original annotations for {variant_record}"
        )
        return variants

    def _expand_annotations_for_all_drugs(
        self,
        original_annotations: List['Annotation'],
        gene_associated_drugs: List[str],
        gene_id: str
    ) -> set['Annotation']:
        """Expands annotations to include all relevant drug associations.

        Creates synthetic annotations for:
        1. Gene-associated drugs not in original annotations
        2. All known drugs from GeneDatabase

        Args:
            original_annotations: Original annotations from TBProfiler
            gene_associated_drugs: Drugs associated with this variant's gene
            gene_id: Gene identifier for database lookup
        Returns:
            Set of all annotations (original + synthetic)
        """
        annotation_set = set(original_annotations)
        seen_drugs = {anno.drug for anno in original_annotations}

        # Add gene-associated drugs
        annotation_set.update(
            self._create_synthetic_annotations(
                base_annotations=original_annotations,
                new_drugs=set(gene_associated_drugs) - seen_drugs,
                confidence="No WHO annotation",
                comment=""
            )
        )
        seen_drugs.update(gene_associated_drugs)

        # Add all drugs from GeneDatabase
        annotation_set.update(
            self._create_synthetic_annotations(
                base_annotations=original_annotations,
                new_drugs=set(GeneDatabase.get_instance().get_drugs(gene_id)) - seen_drugs,
                confidence="No WHO annotation",
                source="Mutation effect for given drug is not in TBDB",
                comment=""
            )
        )

        return annotation_set

    def _create_synthetic_annotations(
        self,
        base_annotations: List['Annotation'],
        new_drugs: set[str],
        confidence: Optional[str] = None,
        source: Optional[str] = None,
        comment: Optional[str] = None,
    ) -> List['Annotation']:
        """Creates synthetic annotations for new drugs based on existing annotation patterns.

        Args:
            base_annotations: Existing annotations to use as templates
            new_drugs: Drugs to create annotations for
            confidence: Confidence level for synthetic annotations
            source: Source information for synthetic annotations
            comment: Comment for synthetic annotations (None = use template's value)
        Returns:
            List of synthetic annotations
        """
        synthetic = []
        for drug in new_drugs:
            for template in base_annotations:
                # Only override non-None values
                overrides = {
                    'drug': drug,
                    'confidence': confidence,
                    'source': source,
                    'comment': comment
                }
                updates = {k: v for k, v in overrides.items() if v is not None}
                synthetic.append(template.model_copy(update=updates))
        return synthetic

    def deduplicate_variants(self, variants: List[Variant]) -> List[Variant]:
        """Deduplicate variants, keeping the one with best annotation.

        When multiple identical variants exist (same gene, position, drug),
        keeps the one with the highest WHO confidence ranking.

        Args:
            variants: List of Variant objects (may contain duplicates)
        Returns:
            List of deduplicated Variants
        """
        unique_variants = {}

        logger.debug(f"Deduplicating Variant objects based on `gene_id`, `gene_name`, `type`, `nucleotide_change`, and `drug` attributes")
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
        if len(variants) != len(result):
            logger.debug(
                f"Deduplicated {len(variants)} Variants to {len(result)} unique Variants"
            )
        else:
            logger.debug("No duplicate Variants found during deduplication step")
        return result

    def generate_unreported_variants(
        self,
        variants: List[Variant],
        sample_id: str,
    ) -> List[Variant]:
        """Generate Variant objects for unreported gene/drug associations.

        Args:
            variants: List of Variant objects
            sample_id: The sample ID associated with the variants.
        Returns:
            List of Variants for unreported genes
        """
        all_unreported_variants = []
        unique_variants = set([variant.gene_id for variant in variants])
        unreported_variants = [v for v in GeneDatabase.get_instance().GENE_DATABASE.keys() if v not in unique_variants]
        for var in unreported_variants:
            drug_list = GeneDatabase.get_instance().get_drugs(var)
            for drug in drug_list:
                variant = Variant.from_thin_air(
                    sample_id=sample_id,
                    gene_id=GeneDatabase.get_instance().get_locus_tag(var),  # type: ignore
                    gene_name=GeneDatabase.get_instance().get_gene_name(var), # type: ignore
                    drug=drug
                )
                all_unreported_variants.append(variant)

        logger.debug(f"Generated {len(all_unreported_variants)} Variants from unreported gene-drug combinations")
        return all_unreported_variants