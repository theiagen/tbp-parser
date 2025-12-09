import pandas as pd
import re

# run at the start and then use as a class feature
def get_position(mutation):
    """  
    This function recieves a mutation (e.g. 'p.Met291Ile') and
    returns the position (numerical part) as an Int
    """    
    pattern = r"-?\d+"

    match = re.findall(pattern, mutation)
    if len(match) > 0:
        return [int(x) for x in match]
    return [0]

def is_within_range(position, range_positions):
    """
    This function determines if a position (aa or nucleotide) are within a particular range
    position is a list of either one or two items
    range_positions is a list of the start and end regions of the range. 
      occasionally, range_positions can be a list of lists, in which case the function will check both ranges and return the result
    """
    if isinstance(range_positions[0], list):
        # check if the value is a list of lists; if so, check both lists
        return is_within_range(position, range_positions[0]) or is_within_range(position, range_positions[1])

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

global ANNOTATION_TO_INTERPRETATION 
ANNOTATION_TO_INTERPRETATION = {
    "Assoc w R": {
        "looker": "R", 
        "mdl": "R"
    },
    "Assoc w R - interim": {
        "looker": "R-Interim",
        "mdl": "R"
    },
    "Assoc w R - Interim": {
        "looker": "R-Interim",
        "mdl": "R"
    },
    "Uncertain significance": {
        "looker": "U", 
        "mdl": "U" 
    },
    "Not assoc w R": {
        "looker": "S",
        "mdl": "S"
    }, 
    "Not assoc w R - Interim": {
        "looker": "S-Interim", 
        "mdl": "S"
    }                              
}
"""
A dictionary to turn TBProfiler WHO annotations into their corresponding Looker or MDL
interpretations; MDL interpretations are the same as Looker but drop the "- Interim" 
designations.
"""

global ANTIMICROBIAL_CODE_TO_GENES
ANTIMICROBIAL_CODE_TO_GENES = {
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
        "STR": {
            "bacA": "STR_bacA",
            "gid": "STR_gid",
            "glpK": "STR_glpK",
            "rpsL": "STR_rpsL",
            "rrs": "STR_rrs",
            "Rv1258c": "STR_Rv1258c",
            "Rv2477c": "STR_Rv2477c",
            "whiB7": "STR_whiB7",
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

global AVERAGE_LOCI_COVERAGE
AVERAGE_LOCI_COVERAGE = {}
"""
A dictionary that will contain the average depth of coverage for each gene
"""

global COVERAGE_DICTIONARY
COVERAGE_DICTIONARY = {}
"""
A dictionary that will contain the percent coverage for each gene
"""

global COVERAGE_THRESHOLD
COVERAGE_THRESHOLD = 100
"""
The coverage threshold (the minimum percent coverage of a gene/locus over the minimum depth)
See also MIN_DEPTH
"""

global DF_LABORATORIAN
DF_LABORATORIAN = pd.DataFrame(columns = [
    "sample_id", "tbprofiler_gene_name", "tbprofiler_locus_tag", "tbprofiler_variant_substitution_type", 
    "tbprofiler_variant_substitution_nt", "tbprofiler_variant_substitution_aa", "confidence", "antimicrobial",
    "looker_interpretation", "mdl_interpretation", "depth", "frequency", "read_support", "rationale", "warning"
])
"""
This dataframe is made to host the laboratorian report
"""

global ETHA237_FREQUENCY
ETHA237_FREQUENCY = 0.1
"""
The minimum frequency for a mutation to be considered for ethA at protein position 237
"""

# to keep - MAYBE CONVERT TO A FILE THOUGH
global GENE_TO_ANTIMICROBIAL_DRUG_NAME
GENE_TO_ANTIMICROBIAL_DRUG_NAME = {
    "aftB": ["amikacin", "capreomycin"],
    "ahpC": ["isoniazid"],
    "ald": ["cycloserine"],
    "alr": ["cycloserine"],
    "atpE": ["bedaquiline"],
    "ccsA": ["amikacin", "capreomycin"],
    "clpC1": ["pyrazinamide"],
    "ddn": ["delamanid"],
    "eis": ["amikacin", "kanamycin"],
    "embA": ["ethambutol"],
    "embB": ["ethambutol"],
    "embC": ["ethambutol"],
    "embR": ["ethambutol"],
    "ethA": ["ethionamide"],
    "ethR": ["ethionamide"],
    "fabG1": ["ethionamide", "isoniazid"],
    "fbiA": ["delamanid"],
    "fbiB": ["delamanid"],
    "fbiC": ["delamanid"],
    "Rv2983": ["delamanid"],
    "fgd1": ["delamanid"],
    "folC": ["para-aminosalicylic_acid"],
    "fprA": ["amikacin", "capreomycin"],
    "gid": ["streptomycin"],
    "gyrA": ["levofloxacin", "moxifloxacin"],
    "gyrB": ["levofloxacin", "moxifloxacin"],
    "inhA": ["ethionamide", "isoniazid"],
    "kasA": ["isoniazid"],
    "katG": ["isoniazid"],
    "mmpL5": ["bedaquiline", "clofazimine"],
    "mmpS5": ["bedaquiline", "clofazimine"],
    "mshA": ["ethionamide", "isoniazid"],
    "ndh": ["ethionamide", "isoniazid"],
    "panD": ["pyrazinamide"],
    "pepQ": ["bedaquiline", "clofazimine"],
    "pncA": ["pyrazinamide"],
    "PPE35": ["pyrazinamide"],
    "ribD": ["para-aminosalicylic_acid"],
    "rplC": ["linezolid"],
    "rpoA": ["rifampin"],
    "rpoB": ["rifampin"],
    "rpoC": ["rifampin"],
    "rpsA": ["pyrazinamide"],
    "rpsL": ["streptomycin"],
    "rrl": ["linezolid"],
    "rrs": ["amikacin", "capreomycin", "kanamycin", "linezolid"],
    "Rv0678": ["bedaquiline", "clofazimine"],
    "Rv1258c": ["isoniazid",  "pyrazinamide", "streptomycin"],
    "Rv1979c": ["bedaquiline", "clofazimine"],
    "Rv2752c": ["isoniazid", "rifampin"],
    "Rv3083": ["ethionamide"],
    "Rv3236c": ["pyrazinamide"],
    "thyA": ["para-aminosalicylic_acid"],
    "thyX": ["para-aminosalicylic_acid"],
    "tlyA": ["capreomycin"],
    "ubiA": ["ethambutol"],
    "whiB6": ["amikacin", "capreomycin", "streptomycin"],
    "whiB7": ["amikacin", "kanamycin", "streptomycin"]
}
"""
A dictionary corresponding each gene to the drug they may confer
resistance to, including the genes in the TBDB watchlist.
See also: https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv, and
https://github.com/jodyphelan/tbdb/blob/master/tbdb.watchlist.csv;
used to create the Laboratorian report
"""

# to keep - MAYBE CONVERT TO A FILE THOUGH
global GENE_TO_LOCUS_TAG
GENE_TO_LOCUS_TAG = {
    "aftB": "Rv3805c",
    "ahpC": "Rv2428",
    "ald": "Rv2780",
    "alr": "Rv3423c",
    "atpE": "Rv1305",
    "ccsA": "Rv0529",
    "clpC1": "Rv3596c",
    "ddn": "Rv3547",
    "eis": "Rv2416c",
    "embA": "Rv3794",
    "embB": "Rv3795",
    "embC": "Rv3793",
    "embR": "Rv1267c",
    "ethA": "Rv3854c",
    "ethR": "Rv3855",
    "fabG1": "Rv1483",
    "fbiA": "Rv3261",
    "fbiB": "Rv3262",
    "fbiC": "Rv1173",
    "Rv2983": "Rv2983",
    "fgd1": "Rv0407",
    "folC": "Rv2447c",
    "fprA": "Rv3106",
    "gid": "Rv3919c",
    "gyrA": "Rv0006",
    "gyrB": "Rv0005",
    "inhA": "Rv1484",
    "kasA": "Rv2245",
    "katG": "Rv1908c",
    "mmpL5": "Rv0676c",
    "mmpS5": "Rv0677c",
    "mshA": "Rv0486",
    "ndh": "Rv1854c",
    "panD": "Rv3601c",
    "pepQ": "Rv2535c",
    "pncA": "Rv2043c",
    "PPE35": "Rv1918c",
    "ribD": "Rv2671",
    "rplC": "Rv0701",
    "rpoA": "Rv3457c",
    "rpoB": "Rv0667",
    "rpoC": "Rv0668",
    "rpsA": "Rv1630",
    "rpsL": "Rv0682",
    "rrl": "EBG00000313339",
    "rrs": "EBG00000313325",
    "Rv0678" : "Rv0678",
    "Rv1258c": "Rv1258c",
    "Rv1979c": "Rv1979c",
    "Rv2752c": "Rv2752c",
    "Rv3083": "Rv3083",
    "Rv3236c": "Rv3236c",
    "thyA": "Rv2764c",
    "thyX": "Rv2754c",
    "tlyA": "Rv1694",
    "ubiA": "Rv3806c",
    "whiB6": "Rv3862c",
    "whiB7": "Rv3197A"
}
"""
A dictionary that matches each gene to its corresponding locus tag;
see also: https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb.bed;
used to create the Laboratorian report
"""  

# to keep - MAYBE CONVERT TO A FILE THOUGH
global GENE_TO_TIER
GENE_TO_TIER = {
    "aftB": "Tier 2", "ahpC": "Tier 1", "atpE": "Tier 1", "ccsA": "Tier 2",
    "clpC1": "Tier 1", "ddn": "Tier 1", "eis": "Tier 1", "embA": "Tier 1",
    "embB": "Tier 1", "embC": "Tier 1", "embR": "Tier 2", "ethA": "Tier 1",
    "ethR": "Tier 2", "fabG1": "Tier 1", "fbiA": "Tier 1", "fbiB": "Tier 1",
    "fbiC": "Tier 1", "fgd1": "Tier 1", "fprA": "Tier 2", "gid": "Tier 1",
    "gyrA": "Tier 1", "gyrB": "Tier 1", "inhA": "Tier 1", "katG": "Tier 1",
    "mmpL5": "Tier 1", "mmpS5": "Tier 1", "mshA": "Tier 2", "ndh": "Tier 2",
    "panD": "Tier 1", "pepQ": "Tier 1", "pncA": "Tier 1", "PPE35": "Tier 2",
    "rplC": "Tier 1", "rpoA": "Tier 2", "rpoB": "Tier 1", "rpoC": "Tier 2",
    "rpsL": "Tier 1", "rrl": "Tier 1", "rrs": "Tier 1", "Rv0678": "Tier 1",
    "Rv1258c": "Tier 1", "Rv1979c": "Tier 2", "Rv2752c": "Tier 2",
    "Rv2983": "Tier 1", "Rv3083": "Tier 2", "Rv3236c": "Tier 2",
    "tlyA": "Tier 1", "ubiA": "Tier 2", "whiB6": "Tier 2", "whiB7": "Tier 1"
}
"""
A dictionary that matches each gene to it's corresponding tier;
tier information provided by CDPH; used to create the Laboratorian report
"""   


global GENES_REPORTED
GENES_REPORTED = set()
"""
A set of genes that have been reported so far
"""

global GENES_WITH_DELETIONS
GENES_WITH_DELETIONS = set()
"""
A set that contains genes that have deletions
"""

global LINEAGE
LINEAGE = ""
"""
A string that contains the lineage of the sample from the "main_lineage" field
"""

global LINEAGE_ENGLISH
LINEAGE_ENGLISH = ""
"""
A string that contains the lineage of the sample in English
"""

global LOW_DEPTH_OF_COVERAGE_LIST
LOW_DEPTH_OF_COVERAGE_LIST = []
"""
A list that will contain the names of genes that have coverages
below the coverage threshold (see also COVERAGE_THRESHOLD)
"""  

global MIN_DEPTH
MIN_DEPTH = 10
"""
The minimum depth threshold
"""

global MIN_FREQUENCY
MIN_FREQUENCY = 0.1
"""
The minimum frequency for a mutation to pass QC
"""

global MIN_LOCUS_PERCENTAGE
MIN_LOCUS_PERCENTAGE = 0.7
"""
The minimum percentage of LIMS genes to pass QC 
for MTBC identification to occur
"""

global MIN_READ_SUPPORT
MIN_READ_SUPPORT = 10
"""
The minimum read support for a mutation to pass QC
(calculated as DEPTH * FREQUENCY)
"""

global MUTATION_FAIL_LIST
MUTATION_FAIL_LIST = []
""" 
A list that will contain the nucleotide mutations that have failed
QC: min coverage < 10x, min freq < 10%, min read support < 10x (or whatever minimums the user specifies)
"""

global OPERATOR
OPERATOR = ""
""" 
This variable holds the operator's name
"""

# to keep - MAYBE CONVERT TO A FILE THOUGH
global WHOV2_PROMOTER_REGIONS
WHOV2_PROMOTER_REGIONS = {
    "aftB": [-129, -1],
    "ahpC": [-93, -1],
    "atpE": [-51, -1],
    "bacA": [-81, -1],
    "ccsA": [-191, -1],
    "clpC1": [-106, -1],
    "ddn": [-51, -1],
    "dnaA": [-314, -1],
    "eis": [-84, -1],
    "embA": [-86, -1],
    "embC": [-1982, -1],
    "embR": [-103, -1],
    "ethA": [-51, -1],
    "ethR": [-26, -1],
    "fbiA": [-138, -1],
    "fbiC": [-127, -1],
    "fgd1": [-51, -1],
    "gid": [-79, -1],
    "glpK": [-52, -1],
    "gyrA": [-35, -1],
    "gyrB": [-108, -1],
    "hadA": [-51, -1],
    "inhA": [-813, -1],
    "katG": [-532, -1],
    "mmpS5": [-85, -1],
    "mshA": [-669, -1],
    "mtrA": [-376, -1],
    "mtrB": [-50, -1],
    "ndh": [-96, -1],
    "nusG": [-201, -1],
    "panD": [[-51, -1], [-1949, -1838]],
    "pepQ": [-51, -1],
    "pncA": [-51, -1],
    "PPE35": [-122, -1],
    "rplC": [[-51, -1], [-503, -323]],
    "rpoA": [-536, -1],
    "rpoB": [-263, -1],
    "rpoC": [-45, -1],
    "rpsA": [-100, -1],
    "rpsL": [-234, -1],
    "rrl": [-51, -1],
    "rrs": [-151, -1],
    "Rv0010c": [-156, -1],
    "Rv0565c": [-78, -1],
    "Rv1129c": [-51, -1],
    "Rv1258c": [-58, -1],
    "Rv1979c": [-470, -1],
    "Rv2477c": [-88, -1],
    "Rv2680": [-153, -1],
    "Rv2681": [-2, -1],
    "Rv2752c": [[-51, -1], [-984, -934]],
    "Rv2983": [-51, -1],
    "Rv3083": [-51, -1],
    "Rv3236c": [[-51, -1], [-538, -488]],
    "sigE": [-51, -1],
    "tlyA": [[-51, -1], [-236, -185]],
    "tsnR": [-51, -1],
    "ubiA": [-51, -1],
    "whiB6": [-126, -1],
    "whiB7": [-404, -1]
}
"""
A dictionary that matches the promoter regions from the WHO v2 catalogue.
If a mutation is within these regions, it needs special consideration; these are nucleotide positions
"""

global RPOB449_FREQUENCY
RPOB449_FREQUENCY = 0.1
"""
The minimum frequency for a mutation to be considered for rpoB at protein position 449
"""

global RRL_FREQUENCY
RRL_FREQUENCY = 0.1
"""
The minimum frequency for a mutation to be considered for rrl
"""

global RRL_READ_SUPPORT
RRL_READ_SUPPORT = 10
"""
The minimum read support for a mutation to be considered for rrl
"""

global RRS_FREQUENCY
RRS_FREQUENCY = 0.1
"""
The minimum frequency for a mutation to be considered for rrs
"""            

global RRS_READ_SUPPORT
RRS_READ_SUPPORT = 10
"""
The minimum read support for a mutation to be considered for rrs
"""

global SEQUENCING_METHOD
SEQUENCING_METHOD = ""
"""
This variable holds the sequencing method
"""

global TNGS_REGIONS
TNGS_REGIONS = {}
"""
The specific primer regions for the gene(s) [tNGS only]
Left blank to not affect WGS sequencing
"""

global TNGS_READ_SUPPORT_BOUNDARIES
TNGS_READ_SUPPORT_BOUNDARIES = []
"""
read support boundaries for tNGS QC; [lower_rs, upper_rs]
"""

global TNGS_FREQUENCY_BOUNDARIES
TNGS_FREQUENCY_BOUNDARIES = []
"""
frequency boundaries for tNGS QC; [lower_f, upper_f]
"""

global TNGS
TNGS = False
"""
A flag to indicate if --tngs is activated
"""

global TREAT_R_AS_S
TREAT_R_AS_S = False
"""
A flag to indicate if --treat-r-as-s is activated; indicates that r mutations in loci 
that are below the coverage threshold will be treated like S or U mutations and will 
not be reported regardless of mutation quality
"""