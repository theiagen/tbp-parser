import logging

logger = logging.getLogger(__name__)

class GeneDatabase:
    GENE_DATABASE = {
        "Rv0005": {
            "locus_tag": "Rv0005",
            "gene_name": "gyrB",
            "tier": "Tier 1",
            "promoter_region": [-108, -1],
            "drugs": ["levofloxacin", "moxifloxacin"]
        },
        # fake fake fake
        "Rv1819c": {
            "locus_tag": "Rv1819c",
            "gene_name": "bacA",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0006": {
            "locus_tag": "Rv0006",
            "gene_name": "gyrA",
            "tier": "Tier 1",
            "promoter_region": [-35, -1],
            "drugs": ["levofloxacin", "moxifloxacin"]
        },
        "Rv0407": {
            "locus_tag": "Rv0407",
            "gene_name": "fgd1",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["delamanid"]
        },
        "Rv0486": {
            "locus_tag": "Rv0486",
            "gene_name": "mshA",
            "tier": "Tier 2",
            "promoter_region": [-669, -1],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv0529": {
            "locus_tag": "Rv0529",
            "gene_name": "ccsA",
            "tier": "Tier 2",
            "promoter_region": [-191, -1],
            "drugs": ["amikacin", "capreomycin"]
        },
        "Rv0667": {
            "locus_tag": "Rv0667",
            "gene_name": "rpoB",
            "tier": "Tier 1",
            "promoter_region": [-263, -1],
            "drugs": ["rifampin"]
        },
        "Rv0668": {
            "locus_tag": "Rv0668",
            "gene_name": "rpoC",
            "tier": "Tier 2",
            "promoter_region": [-45, -1],
            "drugs": ["rifampin"]
        },
        "Rv0676c": {
            "locus_tag": "Rv0676c",
            "gene_name": "mmpL5",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0677c": {
            "locus_tag": "Rv0677c",
            "gene_name": "mmpS5",
            "tier": "Tier 1",
            "promoter_region": [-85, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0678": {
            "locus_tag": "Rv0678",
            "gene_name": "mmpR5",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0682": {
            "locus_tag": "Rv0682",
            "gene_name": "rpsL",
            "tier": "Tier 1",
            "promoter_region": [-234, -1],
            "drugs": ["streptomycin"]
        },
        "Rv0701": {
            "locus_tag": "Rv0701",
            "gene_name": "rplC",
            "tier": "Tier 1",
            "promoter_region": [[-51, -1], [-503, -323]],
            "drugs": ["linezolid"]
        },
        "Rv1173": {
            "locus_tag": "Rv1173",
            "gene_name": "fbiC",
            "tier": "Tier 1",
            "promoter_region": [-127, -1],
            "drugs": ["delamanid"]
        },
        "Rv1258c": {
            "locus_tag": "Rv1258c",
            "gene_name": "Rv1258c",
            "tier": "Tier 1",
            "promoter_region": [-58, -1],
            "drugs": ["isoniazid", "pyrazinamide", "streptomycin"]
        },
        "Rv1267c": {
            "locus_tag": "Rv1267c",
            "gene_name": "embR",
            "tier": "Tier 2",
            "promoter_region": [-103, -1],
            "drugs": ["ethambutol"]
        },
        "Rv1305": {
            "locus_tag": "Rv1305",
            "gene_name": "atpE",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["bedaquiline"]
        },
        "EBG00000313325": {
            "locus_tag": "EBG00000313325",
            "gene_name": "rrs",
            "tier": "Tier 1",
            "promoter_region": [-151, -1],
            "drugs": ["amikacin", "capreomycin", "kanamycin", "linezolid"]
        },
        "EBG00000313339": {
            "locus_tag": "EBG00000313339",
            "gene_name": "rrl",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["linezolid"]
        },
        "Rv1483": {
            "locus_tag": "Rv1483",
            "gene_name": "fabG1",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv1484": {
            "locus_tag": "Rv1484",
            "gene_name": "inhA",
            "tier": "Tier 1",
            "promoter_region": [-813, -1],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv1630": {
            "locus_tag": "Rv1630",
            "gene_name": "rpsA",
            "tier": "NA",
            "promoter_region": [-100, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv1694": {
            "locus_tag": "Rv1694",
            "gene_name": "tlyA",
            "tier": "Tier 1",
            "promoter_region": [[-51, -1], [-236, -185]],
            "drugs": ["capreomycin"]
        },
        "Rv1854c": {
            "locus_tag": "Rv1854c",
            "gene_name": "ndh",
            "tier": "Tier 2",
            "promoter_region": [-96, -1],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv1908c": {
            "locus_tag": "Rv1908c",
            "gene_name": "katG",
            "tier": "Tier 1",
            "promoter_region": [-532, -1],
            "drugs": ["isoniazid"]
        },
        "Rv1918c": {
            "locus_tag": "Rv1918c",
            "gene_name": "PPE35",
            "tier": "Tier 2",
            "promoter_region": [-122, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv1979c": {
            "locus_tag": "Rv1979c",
            "gene_name": "Rv1979c",
            "tier": "Tier 2",
            "promoter_region": [-470, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv2043c": {
            "locus_tag": "Rv2043c",
            "gene_name": "pncA",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv2245": {
            "locus_tag": "Rv2245",
            "gene_name": "kasA",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["isoniazid"]
        },
        "Rv2416c": {
            "locus_tag": "Rv2416c",
            "gene_name": "eis",
            "tier": "Tier 1",
            "promoter_region": [-84, -1],
            "drugs": ["amikacin", "kanamycin"]
        },
        "Rv2428": {
            "locus_tag": "Rv2428",
            "gene_name": "ahpC",
            "tier": "Tier 1",
            "promoter_region": [-93, -1],
            "drugs": ["isoniazid"]
        },
        "Rv2447c": {
            "locus_tag": "Rv2447c",
            "gene_name": "folC",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2535c": {
            "locus_tag": "Rv2535c",
            "gene_name": "pepQ",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv2671": {
            "locus_tag": "Rv2671",
            "gene_name": "ribD",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2752c": {
            "locus_tag": "Rv2752c",
            "gene_name": "Rv2752c",
            "tier": "Tier 2",
            "promoter_region": [[-51, -1], [-984, -934]],
            "drugs": ["isoniazid", "rifampin"]
        },
        "Rv2754c": {
            "locus_tag": "Rv2754c",
            "gene_name": "thyX",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2764c": {
            "locus_tag": "Rv2764c",
            "gene_name": "thyA",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2780": {
            "locus_tag": "Rv2780",
            "gene_name": "ald",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["cycloserine"]
        },
        "Rv2983": {
            "locus_tag": "Rv2983",
            "gene_name": "fbiD",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["delamanid"]
        },
        "Rv3083": {
            "locus_tag": "Rv3083",
            "gene_name": "Rv3083",
            "tier": "Tier 2",
            "promoter_region": [-51, -1],
            "drugs": ["ethionamide"]
        },
        "Rv3106": {
            "locus_tag": "Rv3106",
            "gene_name": "fprA",
            "tier": "Tier 2",
            "promoter_region": [],
            "drugs": ["amikacin", "capreomycin"]
        },
        "Rv3197A": {
            "locus_tag": "Rv3197A",
            "gene_name": "whiB7",
            "tier": "Tier 1",
            "promoter_region": [-404, -1],
            "drugs": ["amikacin", "kanamycin", "streptomycin"]
        },
        "Rv3236c": {
            "locus_tag": "Rv3236c",
            "gene_name": "Rv3236c",
            "tier": "Tier 2",
            "promoter_region": [[-51, -1], [-538, -488]],
            "drugs": ["pyrazinamide"]
        },
        "Rv3261": {
            "locus_tag": "Rv3261",
            "gene_name": "fbiA",
            "tier": "Tier 1",
            "promoter_region": [-138, -1],
            "drugs": ["delamanid"]
        },
        "Rv3262": {
            "locus_tag": "Rv3262",
            "gene_name": "fbiB",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["delamanid"]
        },
        "Rv3423c": {
            "locus_tag": "Rv3423c",
            "gene_name": "alr",
            "tier": "NA",
            "promoter_region": [],
            "drugs": ["cycloserine"]
        },
        "Rv3457c": {
            "locus_tag": "Rv3457c",
            "gene_name": "rpoA",
            "tier": "Tier 2",
            "promoter_region": [-536, -1],
            "drugs": ["rifampin"]
        },
        "Rv3547": {
            "locus_tag": "Rv3547",
            "gene_name": "ddn",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["delamanid"]
        },
        "Rv3596c": {
            "locus_tag": "Rv3596c",
            "gene_name": "clpC1",
            "tier": "Tier 1",
            "promoter_region": [-106, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv3601c": {
            "locus_tag": "Rv3601c",
            "gene_name": "panD",
            "tier": "Tier 1",
            "promoter_region": [[-51, -1], [-1949, -1838]],
            "drugs": ["pyrazinamide"]
        },
        "Rv3793": {
            "locus_tag": "Rv3793",
            "gene_name": "embC",
            "tier": "Tier 1",
            "promoter_region": [-1982, -1],
            "drugs": ["ethambutol"]
        },
        "Rv3794": {
            "locus_tag": "Rv3794",
            "gene_name": "embA",
            "tier": "Tier 1",
            "promoter_region": [-86, -1],
            "drugs": ["ethambutol"]
        },
        "Rv3795": {
            "locus_tag": "Rv3795",
            "gene_name": "embB",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["ethambutol"]
        },
        "Rv3805c": {
            "locus_tag": "Rv3805c",
            "gene_name": "aftB",
            "tier": "Tier 2",
            "promoter_region": [-129, -1],
            "drugs": ["amikacin", "capreomycin"]
        },
        "Rv3806c": {
            "locus_tag": "Rv3806c",
            "gene_name": "ubiA",
            "tier": "Tier 2",
            "promoter_region": [-51, -1],
            "drugs": ["ethambutol"]
        },
        "Rv3854c": {
            "locus_tag": "Rv3854c",
            "gene_name": "ethA",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["ethionamide"]
        },
        "Rv3855": {
            "locus_tag": "Rv3855",
            "gene_name": "ethR",
            "tier": "Tier 2",
            "promoter_region": [-26, -1],
            "drugs": ["ethionamide"]
        },
        "Rv3862c": {
            "locus_tag": "Rv3862c",
            "gene_name": "whiB6",
            "tier": "Tier 2",
            "promoter_region": [-126, -1],
            "drugs": ["amikacin", "capreomycin", "streptomycin"]
        },
        "Rv3919c": {
            "locus_tag": "Rv3919c",
            "gene_name": "gid",
            "tier": "Tier 1",
            "promoter_region": [-79, -1],
            "drugs": ["streptomycin"]
        },
    }

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