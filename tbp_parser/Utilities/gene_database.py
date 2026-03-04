import logging
from pathlib import Path
import yaml

logger = logging.getLogger(__name__)

class GeneDatabase:
    # probably should make this so that the database is loaded from an input parameter instead of hard-coding the file
    # also probably should make this so that the database is only loaded once without needing to pass around an instance of the class
    GENE_DATABASE = yaml.safe_load(open(Path(__file__).parent.parent.parent / "data/gene-database_2026-03-03.yml", "r"))
    
    # Build reverse index: gene_name -> locus_tag
    GENE_DATABASE_INVERTED = {data["gene_name"]: data for data in GENE_DATABASE.values()}

    @classmethod
    def _resolve_database(cls, identifier: str) -> dict | None:
        """Resolve an identifier (gene name or locus tag) to its full gene information."""
        if identifier in cls.GENE_DATABASE:
            return cls.GENE_DATABASE[identifier]
        if identifier in cls.GENE_DATABASE_INVERTED:
            return cls.GENE_DATABASE_INVERTED[identifier]
        logger.debug(f"'{identifier}' not found as gene name or locus tag")
        return None

    @classmethod
    def get_gene_name(cls, identifier: str) -> str | None:
        """Get gene name from locus tag or gene name"""
        db = cls._resolve_database(identifier)
        return db["gene_name"] if db else None

    @classmethod
    def get_locus_tag(cls, identifier: str) -> str | None:
        """Get locus tag from gene name or locus tag"""
        db = cls._resolve_database(identifier)
        return db["locus_tag"] if db else None

    @classmethod
    def get_tier(cls, identifier: str) -> str:
        """Get tier from gene name or locus tag"""
        db = cls._resolve_database(identifier)
        return db["tier"] if db else "NA"

    @classmethod
    def get_promoter_region(cls, identifier: str) -> list[int] | list[list[int]]:
        """Get promoter region from gene name or locus tag"""
        db = cls._resolve_database(identifier)
        return db["promoter_region"] if db else []

    @classmethod
    def get_drugs(cls, identifier: str) -> list[str]:
        """Get associated drugs from gene name or locus tag"""
        db = cls._resolve_database(identifier)
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