import logging
from pathlib import Path
import yaml

logger = logging.getLogger(__name__)

class GeneDatabase:
    _instance = None

    def __init__(self, db_path=None):
        if GeneDatabase._instance is not None:
            raise Exception("GeneDatabase is a singleton. Use GeneDatabase.get_instance().")
        
        if db_path is None:
            db_path = Path(__file__).parent.parent.parent / "data/gene-database_2026-03-03.yml"
        self.GENE_DATABASE = yaml.safe_load(open(db_path, "r"))
        self.GENE_DATABASE_INVERTED = {data["gene_name"]: data for data in self.GENE_DATABASE.values()}
        GeneDatabase._instance = self

    @classmethod
    def get_instance(cls, db_path=None):
        if cls._instance is None:
            cls(db_path)
        return cls._instance

    def _resolve_database(self, identifier: str) -> dict | None:
        """Resolve an identifier (gene name or locus tag) to its full gene information."""
        if identifier in self.GENE_DATABASE:
            return self.GENE_DATABASE[identifier]
        if identifier in self.GENE_DATABASE_INVERTED:
            return self.GENE_DATABASE_INVERTED[identifier]
        logger.debug(f"'{identifier}' not found as gene name or locus tag")
        return None

    def get_gene_name(self, identifier: str) -> str | None:
        """Get gene name from locus tag or gene name"""
        db = self._resolve_database(identifier)
        return db["gene_name"] if db else None

    def get_locus_tag(self, identifier: str) -> str | None:
        """Get locus tag from gene name or locus tag"""
        db = self._resolve_database(identifier)
        return db["locus_tag"] if db else None

    def get_tier(self, identifier: str) -> str:
        """Get tier from gene name or locus tag"""
        db = self._resolve_database(identifier)
        return db["tier"] if db else "NA"

    def get_promoter_region(self, identifier: str) -> list[int] | list[list[int]]:
        """Get promoter region from gene name or locus tag"""
        db = self._resolve_database(identifier)
        return db["promoter_region"] if db else []

    def get_drugs(self, identifier: str) -> list[str]:
        """Get associated drugs from gene name or locus tag"""
        db = self._resolve_database(identifier)
        return db["drugs"] if db else []

    @staticmethod
    def get_gene_to_drug_map(lims_format: dict) -> dict[str, list[str]]:
        """Build a reverse mapping from gene_name to list of drug names from
        a LIMS report format dictionary.

        Args:
            lims_format: dict loaded from LIMS report format YAML
                         (drug -> {code -> {gene -> column}})

        Returns:
            dict[str, list[str]]: gene_name -> list of drug names
        """
        gene_to_drugs = {}
        for drug_name, drug_info in lims_format.items():
            for _code, genes in drug_info.items():
                for gene_name in genes:
                    if gene_name not in gene_to_drugs:
                        gene_to_drugs[gene_name] = []
                    if drug_name not in gene_to_drugs[gene_name]:
                        gene_to_drugs[gene_name].append(drug_name)
        return gene_to_drugs