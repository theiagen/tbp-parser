import json
import logging
from typing import List
from copy import deepcopy

from models.variant import Variant
from utils.gene_database import GeneDatabase

logger = logging.getLogger(__name__)

def json_entry_str(json_entry):
    return f"[{json_entry['gene_name']}][{json_entry['pos']}][{json_entry['type']}][{json_entry['nucleotide_change']}]"

def parse_tbprofiler_json(json_path: str) -> List[Variant]:
    """Parse TBProfiler JSON file into Variant objects.

    Reads the JSON, expands consequences for each variant entry,
    and extracts individual Variant objects (one per drug annotation).

    Args:
        json_path: Path to TBProfiler JSON output file

    Returns:
        List of Variant objects extracted from the JSON
    """
    logger.debug(f"Parsing TBProfiler JSON: {json_path}")

    with open(json_path) as f:
        data = json.load(f)

    logger.debug(f"Parsing `dr_variants` and `other_variants` from JSON")
    logger.debug(f"Found {len(data.get('dr_variants', []))} dr_variants and {len(data.get('other_variants', []))} other_variants")

    # process both dr_variants and other_variants
    all_variant_records = []
    for variant_record in data.get("dr_variants", []) + data.get("other_variants", []):
        all_variant_records += _expand_consequences(variant_record)

    logger.debug(f"Expanded to {len(all_variant_records)} total variant records after parsing consequences")

    # extract variant annotations from each expanded entry
    all_variants = []
    for variant_record in all_variant_records:
        variants = _extract_variant_annotations(variant_record)
        all_variants.extend(variants)

    logger.debug(f"Total Variant objects from JSON: {len(all_variants)}")
    return all_variants


def _expand_consequences(variant_record: dict) -> List[dict]:
    """Expand a variant entry into multiple entries based on consequences.

    Some variants affect multiple genes (e.g., a mutation in mmpR5 that
    also affects mmpL5 and mmpS5). This creates separate entries for each
    affected gene.

    Args:
        variant_record: Single variant entry from TBProfiler JSON

    Returns:
        List of variant entries (original + consequences)
    """
    # Only expand if consequences exist and gene is one of mmpR5/mmpL5/mmpS5 otherwise return original
    if (
        not variant_record.get("consequences") or
        variant_record.get('gene_id') not in ['Rv0676c', 'Rv0677c', 'Rv0678']
    ):
        return [variant_record]

    logger.debug(f"Expanding {len(variant_record['consequences'])} consequences for JSON entry: {json_entry_str(variant_record)}")

    all_variant_records = []
    for consequence in variant_record['consequences']:
        # Create a copy and update with consequence-specific data
        new_entry = deepcopy(variant_record)
        new_entry.update({
            'gene_id': consequence['gene_id'],
            'gene_name': consequence['gene_name'],
            'feature_id': consequence['feature_id'],
            'type': consequence['type'],
            'nucleotide_change': consequence['nucleotide_change'],
            'protein_change': consequence['protein_change'],
            'annotation': consequence['annotation'],
        })
        all_variant_records.append(new_entry)

    #logger.debug(f"Consequence entries: {[json_entry_str(entry) for entry in all_variant_records]}")

    return all_variant_records


def _extract_variant_annotations(variant_record: dict) -> List[Variant]:
    """Extract Variant objects from a single variant entry.

    Creates one Variant per drug annotation, plus additional Variants
    for drugs associated with the gene that weren't explicitly annotated.

    Args:
        variant_record: Single variant entry (possibly expanded from consequences)

    Returns:
        List of Variant objects
    """
    variants_list = []
    seen_drugs = set()

    logger.debug(f"Extracting annotations for JSON entry: {json_entry_str(variant_record)}")

    for annotation in variant_record.get("annotation", []):
        # normalize confidence, maybe this should live somewhere else?
        confidence = annotation['confidence']
        if annotation.get('comment') == "Not found in WHO catalogue":
            confidence = "No WHO annotation"

        variant = Variant(
            pos=variant_record['pos'],
            depth=variant_record['depth'],
            freq=variant_record['freq'],
            gene_id=variant_record['gene_id'],
            gene_name=variant_record['gene_name'],
            type=variant_record['type'],
            nucleotide_change=variant_record['nucleotide_change'],
            protein_change=variant_record['protein_change'],
            confidence=confidence,
            drug=annotation['drug'],
            source=annotation['source'],
            comment=annotation['comment'],
        )
        variants_list.append(variant)
        seen_drugs.add(annotation['drug'])

    # add variants for drugs associated with the gene but not in annotations
    gene_associated_drug_list = set(variant_record.get("gene_associated_drugs", [])) - seen_drugs
    for drug in gene_associated_drug_list:
        source = "Mutation effect for given drug is not in TBDB"
        variants_list.append(Variant(
            pos=variant_record['pos'],
            depth=variant_record['depth'],
            freq=variant_record['freq'],
            gene_id=variant_record['gene_id'],
            gene_name=variant_record['gene_name'],
            type=variant_record['type'],
            nucleotide_change=variant_record['nucleotide_change'],
            protein_change=variant_record['protein_change'],
            confidence="No WHO annotation",
            drug=drug,
            source=source,
            comment="",
        ))
        seen_drugs.add(drug)

    # add variants for any drugs in the gene database not yet seen
    remaining_drug_list = set(GeneDatabase.get_drugs(variant_record['gene_id'])) - seen_drugs
    for drug in remaining_drug_list:
        variants_list.append(Variant(
            pos=variant_record['pos'],
            depth=variant_record['depth'],
            freq=variant_record['freq'],
            gene_id=variant_record['gene_id'],
            gene_name=variant_record['gene_name'],
            type=variant_record['type'],
            nucleotide_change=variant_record['nucleotide_change'],
            protein_change=variant_record['protein_change'],
            confidence="No WHO annotation",
            drug=drug,
            source="",
            comment="",
        ))
        seen_drugs.add(drug)

    # this is all just for logging/debugging
    variant_sources = []
    annotation_drugs = [_.get('drug') for _ in variant_record.get('annotation', [])]
    if annotation_drugs:
        variant_sources.append(f"`annotations`: {annotation_drugs}")
    if gene_associated_drug_list:
        variant_sources.append(f"`gene_associated_drugs`: {list(gene_associated_drug_list)}")
    if remaining_drug_list:
        variant_sources.append(f"`GeneDatabase`: {list(remaining_drug_list)}")

    logger.debug(f"Adding {len(variants_list)} Variant objects from: {', '.join(variant_sources) if variant_sources else 'N/A'}")
    return variants_list
