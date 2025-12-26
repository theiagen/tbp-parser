import re

def get_position(mutation) -> list[int]:
    """This function recieves a mutation and returns the position as an integer

    Args:
        mutation (str): the mutation; can be either p.Met291Ile or p.Lys123_delArg125

    Returns:
        list[int]: the numerical position(s) of the mutation (in the above example, [291] or [123, 125])
    """    
    pattern = r"-?\d+"
    match = re.findall(pattern, mutation)
    if len(match) > 0:
        return [int(x) for x in match]
    return [None]

def get_mutation_genomic_positions(position, mutation) -> tuple[int, int]:
    """This function receives the genomic position and a mutation and returns the genomic position range as a list of integers

    Args:
        position(int): the genomic position of the mutation (e.g., 2000)
        mutation (str): the nucleotide mutation; can be either c.-33_327del or c.1693G>T

    Returns:
        tuple[int, int]: the start and stop genomic positions of the mutation (in the above example, [2000, 2360] or [2000, 2000])
    """    
    pattern = r"-?\d+"
    match = re.findall(pattern, mutation)
    if len(match) == 1:
        return (position, position) 
    elif len(match) == 2:
        return (position, position + (abs(int(match[0]) - int(match[1]))))
    return (None, None)

def is_mutation_within_range(position, range_positions) -> bool:
    """Determines if a position is within a particular range

    Args:
        position (list[int]): either one or two positions
        range_positions (list[int] or list[list[int]]): the start and end regions of the range

    Returns:
        bool: true if the position is within the range_positions, false otherwise
    """    
    try:
        if isinstance(range_positions[0], list):
            # check if the value is a list of lists; if so, check both lists
            return is_mutation_within_range(position, range_positions[0]) or is_mutation_within_range(position, range_positions[1])

        if len(position) > 1:
            # if the value is a list of two items, check if the position is within the range
            if any([x in range(range_positions[0], range_positions[1]) for x in position]):
                return True
            if any([x in range(position[0], position[1]) for x in range_positions]):
                return True

        # the position is a single item
        elif range_positions[0] <= position[0] <= range_positions[1]:
            return True

        return False
    except:
        return False

global DRUG_COLUMNS_TO_GENE_COLUMNS
DRUG_COLUMNS_TO_GENE_COLUMNS = {
    "amikacin": {
        "AMK": {
            "bacA": "AMK_bacA",
            "ccsA": "AMK_ccsA",
            "eis": "AMK_eis",
            "rrs": "AMK_rrs",
            "Rv2477C": "AMK_Rv2477C",
            "whiB6": "AMK_whiB6",
            "whiB7": "AMK_whiB7",
        }
    },
    "bedaquiline": {
        "BDQ": {
            "atpE": "BDQ_atpE",
            "lpqB": "BDQ_lpqB",
            "mmpL5": "BDQ_mmpL5",
            "mmpR5": "BDQ_mmpR5",
            "mmpS5": "BDQ_mmpS5",
            "mtrA": "BDQ_mtrA",
            "mtrB": "BDQ_mtrB",
            "pepQ": "BDQ_pepQ",
            "Rv1979c": "BDQ_Rv1979c",
        }
    },
    "capreomycin": {
        "CAP": {
            "bacA": "CAP_bacA",
            "ccsA": "CAP_ccsA",
            "rrl": "CAP_rrl",
            "rrs": "CAP_rrs",
            "Rv2680": "CAP_Rv2680",
            "Rv2681": "CAP_Rv2681",
            "tlyA": "CAP_tlyA",
            "whiB6": "CAP_whiB6",
        }
    },
    "clofazimine": {
      "CFZ": {
            "fbiA": "CFZ_fbiA",
            "fbiB": "CFZ_fbiB",
            "fbiC": "CFZ_fbiC",
            "fbiD": "CFZ_fbiD",
            "fgd1": "CFZ_fgd1",
            "mmpL5": "CFZ_mmpL5",
            "mmpR5": "CFZ_mmpR5",
            "mmpS5": "CFZ_mmpS5",
            "pepQ": "CFZ_pepQ",
            "Rv1979c": "CFZ_Rv1979c",
        }
    },
    "cycloserine": {
        "CS": {
            "ald": "CS_ald",
            "alr": "CS_alr"
        }
    },
    "delamanid": {
        "DLM": {
            "ddn": "DLM_ddn",
            "fbiA": "DLM_fbiA",
            "fbiB": "DLM_fbiB",
            "fbiC": "DLM_fbiC",
            "fbiD": "DLM_fbiD",
            "fgd1": "DLM_fgd1",
            "ndh": "DLM_ndh",
        }
    },
    "ethambutol": {
        "EMB": {
            "aftB": "EMB_aftB",
            "embA": "EMB_embA",
            "embB": "EMB_embB",
            "embC": "EMB_embC",
            "embR": "EMB_embR",
            "glpK": "EMB_glpK",
            "Rv2477c": "EMB_Rv2477c",
            "Rv2752c": "EMB_Rv2752c",
            "ubiA": "EMB_ubiA",
        }
    },
    "ethionamide": {
        "ETO": {
            "ethA": "ETO_ethA",
            "ethR": "ETO_ethR",
            "inhA": "ETO_inhA",
            "mshA": "ETO_mshA",
            "ndh": "ETO_ndh",
            "Rv0565c": "ETO_Rv0565c",
            "Rv3083": "ETO_Rv3083",
        }
    },
    "isoniazid": {
        "INH": { 
            "ahpC": "INH_ahpC",
            "dnaA": "INH_dnaA",
            "glpK": "INH_glpK",
            "hadA": "INH_hadA",
            "inhA": "INH_inhA",
            "kasA": "INH_kasA",
            "katG": "INH_katG",
            "mmaA3": "INH_mmaA3",
            "mshA": "INH_mshA",
            "ndh": "INH_ndh",
            "Rv0010c": "INH_Rv0010c",
            "Rv1129c": "INH_Rv1129c",
            "Rv1258c": "INH_Rv1258c",
            "Rv2752c": "INH_Rv2752c",
        }
    },
    "kanamycin": {
        "KAN": {
            "bacA": "KAN_bacA",
            "ccsA": "KAN_ccsA",
            "eis": "KAN_eis",
            "rrs": "KAN_rrs",
            "Rv2477c": "KAN_Rv2477c",
            "whiB6": "KAN_whiB6",
            "whiB7": "KAN_whiB7",
        }
    },
    "levofloxacin": {
        "LFX": {
            "glpK": "LFX_glpK",
            "gyrA": "LFX_gyrA",
            "gyrB": "LFX_gyrB",
            "Rv1129c": "LFX_Rv1129c",
            "Rv2477c": "LFX_Rv2477c",
            "Rv2752c": "LFX_Rv2752c",
        }
    },
    "linezolid": {
        "LZD": {
            "rplC": "LZD_rplC",
            "rrl": "LZD_rrl",
            "tsnR": "LZD_tsnR",
        }
    },
    "moxifloxacin": {
        "MFX": {
            "glpK": "MFX_glpK",
            "gyrA": "MFX_gyrA",
            "gyrB": "MFX_gyrB",
            "Rv1129c": "MFX_Rv1129c",
            "Rv2477c": "MFX_Rv2477c",
            "Rv2752c": "MFX_Rv2752c",
        }
    },
    "para-aminosalicylic_acid": {
        "PAS": {
            "folC": "PAS_folC",
            "ribD": "PAS_ribD",
            "thyA": "PAS_thyA",
            "thyX": "PAS_thyX",
        }
    },
    "pretomanid": {
        "PMD": {
            "ddn": "PMD_ddn",
            "fbiA": "PMD_fbiA",
            "fbiB": "PMD_fbiB",
            "fbiC": "PMD_fbiC",
            "fbiD": "PMD_fbiD",
            "fgd1": "PMD_fgd1",
        }
    },
    "pyrazinamide": {
        "PZA": {
            "clpC1": "PZA_clpC1",
            "panD": "PZA_panD",
            "pncA": "PZA_pncA",
            "PPE35": "PZA_PPE35",
            "rpsA": "PZA_rpsA",
            "Rv1258c": "PZA_Rv1258c",
            "Rv3236c": "PZA_Rv3236c",
            "sigE": "PZA_sigE",
        }
    },
    "rifampicin": {
        "RIF": {
            "glpK": "RIF_glpK",
            "lpqB": "RIF_lpqB",
            "mtrA": "RIF_mtrA",
            "mtrB": "RIF_mtrB",
            "nusG": "RIF_nusG",
            "rpoA": "RIF_rpoA",
            "rpoB": "RIF_rpoB",
            "rpoC": "RIF_rpoC",
            "Rv1129c": "RIF_Rv1129c",
            "Rv2477c": "RIF_Rv2477c",
            "Rv2752c": "RIF_Rv2752c",
        }
    },
    "streptomycin": {
        "STM": {
            "bacA": "STM_bacA",
            "gid": "STM_gid",
            "glpK": "STM_glpK",
            "rpsL": "STM_rpsL",
            "rrs": "STM_rrs",
            "Rv1258c": "STM_Rv1258c",
            "Rv2477c": "STM_Rv2477c",
            "whiB7": "STM_whiB7",
        }
    }
}
"""
This relatively complex dictionary that takes the following format:

{
    "antimicrobial_drug_name": {
        "antimicrobial_column_name_in_lims_report": {
            "gene_name": "antimicrobial_gene_combo_column_name_in_lims_report",
            ...
        }
    },
    ...
}

Used to create the LIMS report, this dictionary can be easily extended to include
new drugs/genes or reduced to only the drugs/genes needed in the configuration file.

This is derived from the TBProfiler `genes.bed` file, which can be accessed here: 
https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb/genes.bed.
"""

global OUTPUT_RENAMING
OUTPUT_RENAMING = {}
"""This dictionary maps existing default language to new language for output reports.
For example, if the desired output language for 'rifampicin' is 'rifampin', the dictionary
would contain the entry {'rifampicin': 'rifampin'}. At the end of processing, all instances
of 'rifampicin' in the output reports would be replaced with 'rifampin'.
"""
