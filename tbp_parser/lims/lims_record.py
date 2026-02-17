from typing import Any, Optional
from pydantic import BaseModel, Field
from variant import Variant
from utils import Helper

class LIMSGeneCode(BaseModel):
    """Class representing gene-specific results for the LIMS report."""
    gene_code: str

    # Fields for storing the results populated in lims_processor.py
    gene_target_value: Optional[str] = Field(default=None, exclude=True)
    max_mdl_interpretation: Optional[str] = Field(default=None, exclude=True)
    max_mdl_variants: list[Variant] = Field(default_factory=list, exclude=True)

class LIMSRecord(BaseModel):
    """Class representing drug-specific results for the LIMS report."""
    drug: str
    drug_code: str
    gene_codes: dict[str, LIMSGeneCode]

    # Fields for storing the results populated in lims_processor.py
    drug_target_value: Optional[str] = Field(default=None, exclude=True)

    # Post-init processing to compute derived attributes
    def model_post_init(self, __context: Any = None):
        Helper.normalize_field_values(self)
