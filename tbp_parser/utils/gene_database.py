import logging


class GeneDatabase:
    GENE_DATABASE = {
        "Rv0005": {
            "gene_name": "gyrB",
            "tier": "Tier 1",
            "promoter_region": [-108, -1],
            "drugs": ["levofloxacin", "moxifloxacin"]
        },
        "Rv0006": {
            "gene_name": "gyrA",
            "tier": "Tier 1",
            "promoter_region": [-35, -1],
            "drugs": ["levofloxacin", "moxifloxacin"]
        },
        "Rv0407": {
            "gene_name": "fgd1",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["delamanid"]
        },
        "Rv0486": {
            "gene_name": "mshA",
            "tier": "Tier 2",
            "promoter_region": [-669, -1],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv0529": {
            "gene_name": "ccsA",
            "tier": "Tier 2",
            "promoter_region": [-191, -1],
            "drugs": ["amikacin", "capreomycin"]
        },
        "Rv0667": {
            "gene_name": "rpoB",
            "tier": "Tier 1",
            "promoter_region": [-263, -1],
            "drugs": ["rifampin"]
        },
        "Rv0668": {
            "gene_name": "rpoC",
            "tier": "Tier 2",
            "promoter_region": [-45, -1],
            "drugs": ["rifampin"]
        },
        "Rv0676c": {
            "gene_name": "mmpL5",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0677c": {
            "gene_name": "mmpS5",
            "tier": "Tier 1",
            "promoter_region": [-85, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0678": {
            "gene_name": "mmpR5",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv0682": {
            "gene_name": "rpsL",
            "tier": "Tier 1",
            "promoter_region": [-234, -1],
            "drugs": ["streptomycin"]
        },
        "Rv0701": {
            "gene_name": "rplC",
            "tier": "Tier 1",
            "promoter_region": [[-51, -1], [-503, -323]],
            "drugs": ["linezolid"]
        },
        "Rv1173": {
            "gene_name": "fbiC",
            "tier": "Tier 1",
            "promoter_region": [-127, -1],
            "drugs": ["delamanid"]
        },
        "Rv1258c": {
            "gene_name": "Rv1258c",
            "tier": "Tier 1",
            "promoter_region": [-58, -1],
            "drugs": ["isoniazid", "pyrazinamide", "streptomycin"]
        },
        "Rv1267c": {
            "gene_name": "embR",
            "tier": "Tier 2",
            "promoter_region": [-103, -1],
            "drugs": ["ethambutol"]
        },
        "Rv1305": {
            "gene_name": "atpE",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["bedaquiline"]
        },
        "EBG00000313325": {
            "gene_name": "rrs",
            "tier": "Tier 1",
            "promoter_region": [-151, -1],
            "drugs": ["amikacin", "capreomycin", "kanamycin", "linezolid"]
        },
        "EBG00000313339": {
            "gene_name": "rrl",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["linezolid"]
        },
        "Rv1483": {
            "gene_name": "fabG1",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv1484": {
            "gene_name": "inhA",
            "tier": "Tier 1",
            "promoter_region": [-813, -1],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv1630": {
            "gene_name": "rpsA",
            "tier": "",
            "promoter_region": [-100, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv1694": {
            "gene_name": "tlyA",
            "tier": "Tier 1",
            "promoter_region": [[-51, -1], [-236, -185]],
            "drugs": ["capreomycin"]
        },
        "Rv1854c": {
            "gene_name": "ndh",
            "tier": "Tier 2",
            "promoter_region": [-96, -1],
            "drugs": ["ethionamide", "isoniazid"]
        },
        "Rv1908c": {
            "gene_name": "katG",
            "tier": "Tier 1",
            "promoter_region": [-532, -1],
            "drugs": ["isoniazid"]
        },
        "Rv1918c": {
            "gene_name": "PPE35",
            "tier": "Tier 2",
            "promoter_region": [-122, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv1979c": {
            "gene_name": "Rv1979c",
            "tier": "Tier 2",
            "promoter_region": [-470, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv2043c": {
            "gene_name": "pncA",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv2245": {
            "gene_name": "kasA",
            "tier": "",
            "promoter_region": [],
            "drugs": ["isoniazid"]
        },
        "Rv2416c": {
            "gene_name": "eis",
            "tier": "Tier 1",
            "promoter_region": [-84, -1],
            "drugs": ["amikacin", "kanamycin"]
        },
        "Rv2428": {
            "gene_name": "ahpC",
            "tier": "Tier 1",
            "promoter_region": [-93, -1],
            "drugs": ["isoniazid"]
        },
        "Rv2447c": {
            "gene_name": "folC",
            "tier": "",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2535c": {
            "gene_name": "pepQ",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["bedaquiline", "clofazimine"]
        },
        "Rv2671": {
            "gene_name": "ribD",
            "tier": "",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2752c": {
            "gene_name": "Rv2752c",
            "tier": "Tier 2",
            "promoter_region": [[-51, -1], [-984, -934]],
            "drugs": ["isoniazid", "rifampin"]
        },
        "Rv2754c": {
            "gene_name": "thyX",
            "tier": "",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2764c": {
            "gene_name": "thyA",
            "tier": "",
            "promoter_region": [],
            "drugs": ["para-aminosalicylic_acid"]
        },
        "Rv2780": {
            "gene_name": "ald",
            "tier": "",
            "promoter_region": [],
            "drugs": ["cycloserine"]
        },
        "Rv2983": {
            "gene_name": "fbiD",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["delamanid"]
        },
        "Rv3083": {
            "gene_name": "Rv3083",
            "tier": "Tier 2",
            "promoter_region": [-51, -1],
            "drugs": ["ethionamide"]
        },
        "Rv3106": {
            "gene_name": "fprA",
            "tier": "Tier 2",
            "promoter_region": [],
            "drugs": ["amikacin", "capreomycin"]
        },
        "Rv3197A": {
            "gene_name": "whiB7",
            "tier": "Tier 1",
            "promoter_region": [-404, -1],
            "drugs": ["amikacin", "kanamycin", "streptomycin"]
        },
        "Rv3236c": {
            "gene_name": "Rv3236c",
            "tier": "Tier 2",
            "promoter_region": [[-51, -1], [-538, -488]],
            "drugs": ["pyrazinamide"]
        },
        "Rv3261": {
            "gene_name": "fbiA",
            "tier": "Tier 1",
            "promoter_region": [-138, -1],
            "drugs": ["delamanid"]
        },
        "Rv3262": {
            "gene_name": "fbiB",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["delamanid"]
        },
        "Rv3423c": {
            "gene_name": "alr",
            "tier": "",
            "promoter_region": [],
            "drugs": ["cycloserine"]
        },
        "Rv3457c": {
            "gene_name": "rpoA",
            "tier": "Tier 2",
            "promoter_region": [-536, -1],
            "drugs": ["rifampin"]
        },
        "Rv3547": {
            "gene_name": "ddn",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["delamanid"]
        },
        "Rv3596c": {
            "gene_name": "clpC1",
            "tier": "Tier 1",
            "promoter_region": [-106, -1],
            "drugs": ["pyrazinamide"]
        },
        "Rv3601c": {
            "gene_name": "panD",
            "tier": "Tier 1",
            "promoter_region": [[-51, -1], [-1949, -1838]],
            "drugs": ["pyrazinamide"]
        },
        "Rv3793": {
            "gene_name": "embC",
            "tier": "Tier 1",
            "promoter_region": [-1982, -1],
            "drugs": ["ethambutol"]
        },
        "Rv3794": {
            "gene_name": "embA",
            "tier": "Tier 1",
            "promoter_region": [-86, -1],
            "drugs": ["ethambutol"]
        },
        "Rv3795": {
            "gene_name": "embB",
            "tier": "Tier 1",
            "promoter_region": [],
            "drugs": ["ethambutol"]
        },
        "Rv3805c": {
            "gene_name": "aftB",
            "tier": "Tier 2",
            "promoter_region": [-129, -1],
            "drugs": ["amikacin", "capreomycin"]
        },
        "Rv3806c": {
            "gene_name": "ubiA",
            "tier": "Tier 2",
            "promoter_region": [-51, -1],
            "drugs": ["ethambutol"]
        },
        "Rv3854c": {
            "gene_name": "ethA",
            "tier": "Tier 1",
            "promoter_region": [-51, -1],
            "drugs": ["ethionamide"]
        },
        "Rv3855": {
            "gene_name": "ethR",
            "tier": "Tier 2",
            "promoter_region": [-26, -1],
            "drugs": ["ethionamide"]
        },
        "Rv3862c": {
            "gene_name": "whiB6",
            "tier": "Tier 2",
            "promoter_region": [-126, -1],
            "drugs": ["amikacin", "capreomycin", "streptomycin"]
        },
        "Rv3919c": {
            "gene_name": "gid",
            "tier": "Tier 1",
            "promoter_region": [-79, -1],
            "drugs": ["streptomycin"]
        },
    }

    def __init__(self, locus_tag: str, tier: str, promoter_region: list[int] | list[list[int]], drugs: list[str]) -> None:
        self.logger = logging.getLogger(self.__class__.__name__)

        self.locus_tag = locus_tag
        self.tier = tier
        self.promoter_region = promoter_region
        self.drugs = drugs

# Build reverse index: gene_name -> locus_tag
    _GENE_NAME_TO_LOCUS = {info["gene_name"]: locus_tag
                           for locus_tag, info in GENE_DATABASE.items()}

    @classmethod
    def _resolve_to_locus_tag(cls, identifier: str) -> str:
        """
        Convert gene name or locus tag to locus tag.
        If identifier is already a locus tag, return it.
        If identifier is a gene name, look up its locus tag.
        """
        # Check if it's already a locus tag
        if identifier in cls.GENE_DATABASE:
            return identifier
        # Check if it's a gene name
        if identifier in cls._GENE_NAME_TO_LOCUS:
            return cls._GENE_NAME_TO_LOCUS[identifier]
        # Not found
        raise KeyError(f"'{identifier}' not found as gene name or locus tag")

    @classmethod
    def get_gene_name(cls, identifier: str) -> str:
        """Get gene name from locus tag or gene name"""
        locus_tag = cls._resolve_to_locus_tag(identifier)
        return cls.GENE_DATABASE[locus_tag]["gene_name"]

    @classmethod
    def get_locus_tag(cls, identifier: str) -> str:
        """Get locus tag from gene name or locus tag"""
        return cls._resolve_to_locus_tag(identifier)

    @classmethod
    def get_tier(cls, identifier: str) -> str | None:
        """Get tier from gene name or locus tag"""
        locus_tag = cls._resolve_to_locus_tag(identifier)
        return cls.GENE_DATABASE[locus_tag]["tier"]

    @classmethod
    def get_promoter_region(cls, identifier: str):
        """Get promoter region from gene name or locus tag"""
        locus_tag = cls._resolve_to_locus_tag(identifier)
        return cls.GENE_DATABASE[locus_tag]["promoter_region"]

    @classmethod
    def get_drugs(cls, identifier: str) -> list[str]:
        """Get associated drugs from gene name or locus tag"""
        locus_tag = cls._resolve_to_locus_tag(identifier)
        return cls.GENE_DATABASE[locus_tag]["drugs"]