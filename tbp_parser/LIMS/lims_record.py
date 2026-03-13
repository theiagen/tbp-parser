from typing import Optional
from pydantic import BaseModel, Field
from tbp_parser.Variant.variant import Variant

class LIMSGeneCode(BaseModel):
    """Class representing gene-specific results for the LIMS report."""
    gene_code: str

    # Fields for storing the results populated in lims_processor.py
    gene_target_value: Optional[str] = Field(default=None, exclude=True)
    max_mdl_interpretation: Optional[str] = Field(default=None, exclude=True)
    max_mdl_variants: list[Variant] = Field(default_factory=list, exclude=True)

    def __str__(self):
        return f"LIMSGeneCode(gene_code={self.gene_code}, gene_target_value={self.gene_target_value}, max_mdl_interpretation={self.max_mdl_interpretation}, max_mdl_variants={self.max_mdl_variants})"

    def __repr__(self):
        return f"LIMSGeneCode(gene_code={self.gene_code}, gene_target_value={self.gene_target_value}, max_mdl_interpretation={self.max_mdl_interpretation}, max_mdl_variants={self.max_mdl_variants})"
class LIMSRecord(BaseModel):
    """Class representing drug-specific results for the LIMS report."""
    drug: str
    drug_code: str
    gene_codes: dict[str, LIMSGeneCode]

    # Fields for storing the results populated in lims_processor.py
    drug_target_value: Optional[str] = Field(default=None, exclude=True)

    def __str__(self):
        return f"LIMSRecord(drug={self.drug}, drug_code={self.drug_code}, drug_target_value={self.drug_target_value})"

    def __repr__(self):
        return f"LIMSRecord(drug={self.drug}, drug_code={self.drug_code}, drug_target_value={self.drug_target_value})"