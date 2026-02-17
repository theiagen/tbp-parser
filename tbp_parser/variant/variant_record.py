from typing import Any, List
from utils import Helper
from pydantic import BaseModel
import logging

logger = logging.getLogger(__name__)

class Annotation(BaseModel):
    """Data class representing a single annotation entry (dict) from a list of dicts under the `annotation` field in
    the TBProfiler JSON output. This class is primarily used for organizational purposes to define the structure of an annotation.

    Can be found in the:
    - `annotation` field of a VariantRecord.
    - `annotation` field of a Consequences entry within a VariantRecord.
    """
    drug: str
    confidence: str
    source: str
    comment: str
    model_config = {"extra": "ignore"}

    def __hash__(self):
        # Custom hash method to allow usage in sets and define uniqueness
        return hash((self.drug, self.confidence, self.source, self.comment))

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        Helper.normalize_field_values(self)

class Consequences(BaseModel):
    """Data class representing a single consequence entry (dict) from a list of dicts under the `consequences` field in
    the TBProfiler JSON output. This class is primarily used for organizational purposes to define the structure of a consequence.

    Can be found in the `consequences` field of a VariantRecord. Consequences will often be expanded into multiple VariantRecord entries
    depending on the original gene affected (mmpR5, mmpL5, and mmpS5). Consequences also have their own `annotation` field which may differ from the parent VariantRecord.
    These will be used to create additional Variant objects later on.
    """
    gene_id: str
    gene_name: str
    feature_id: str
    type: str
    nucleotide_change: str
    protein_change: str
    annotation: List[Annotation]
    model_config = {"extra": "ignore"}

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        Helper.normalize_field_values(self)

class VariantRecord(BaseModel):
    """Data class representing a single variant record (aka JSON entry) from TBProfiler JSON output.

    This class is primarily used for internal use to help organize the JSON output. A lot of the fields are never used,
    but are included here to explicitly define the structure of a VariantRecord as represented in the TBProfiler JSON output.

    Each VariantRecord may correspond to and ultimately spawn multiple Variant objects, depending on the:
    - `annotation` field (one Variant per drug annotation)
    - `consequences` field (may expand and create multiple additional VariantRecords each with their own annotations)
    - `gene_associated_drugs` field (may create additional Variants for drugs with no annotations)
    """
    sample_id: str
    pos: int
    depth: int
    freq: float
    gene_id: str
    gene_name: str
    type: str
    nucleotide_change: str
    protein_change: str
    annotation: List[Annotation]
    consequences: List[Consequences]
    gene_associated_drugs: List[str]

    def __str__(self) -> str:
        return f"VariantRecord([{self.gene_name}][{self.pos}][{self.type}][{self.nucleotide_change}])"

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        Helper.normalize_field_values(self)

    @classmethod
    def from_consequences(
        cls,
        variant_record: 'VariantRecord',
        consequences: 'Consequences',
    ) -> 'VariantRecord':
        """Returns a new VariantRecord with the attributes of `consequences` overwriting the corresponding fields of the copied VariantRecord.
        Args:
            variant_record: The original VariantRecord instance.
            consequences: A Consequences instance representing a single consequences entry. (dict not the list of dicts)
        Returns:
            VariantRecord: A copy of the original VariantRecord instance, but populated with the new `consequences` attributes where applicable.
        """
        consequences_dict = consequences.model_dump()
        # preserve original consequence annotations if they exist and are not empty
        # otherwise, create a default annotation based on the gene_associated drugs
        # https://github.com/theiagen/tbp-parser/blob/0dae65710a0f14e0796cb1b7a3a1dc390e2aeb0f/tbp_parser/Variant.py#L64
        if not hasattr(consequences, "annotation") or not consequences.annotation:
            consequences_dict["annotation"] = [Annotation(
                drug=drug,
                confidence="No WHO annotation",
                source="",
                comment="",
            ) for drug in variant_record.gene_associated_drugs]

        return variant_record.model_copy(update=consequences_dict)
