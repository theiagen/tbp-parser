import logging
import yaml

logger = logging.getLogger(__name__)

class GeneDatabase:
    GENE_DATABASE: dict[str, dict] = {}
    GENE_DATABASE_INVERTED: dict[str, dict] = {}

    def __init__(self, db_path: str):
        GeneDatabase.GENE_DATABASE = yaml.safe_load(open(db_path, "r"))
        GeneDatabase.GENE_DATABASE_INVERTED = {data["gene_name"]: data for data in GeneDatabase.GENE_DATABASE.values()}

    @classmethod
    def _resolve_database(cls, identifier: str) -> dict | None:
        """Resolve an identifier (gene name or locus tag) to its full gene information."""
        if identifier in cls.GENE_DATABASE:
            return cls.GENE_DATABASE[identifier]
        if identifier in cls.GENE_DATABASE_INVERTED:
            return cls.GENE_DATABASE_INVERTED[identifier]
        return None

    @classmethod
    def get_db(cls) -> dict:
        """Get the entire gene database."""
        return cls.GENE_DATABASE

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