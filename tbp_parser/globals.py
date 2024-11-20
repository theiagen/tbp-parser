import pandas as pd
import re

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
  print("function")
  print("position", position)
  print("range_positions", range_positions)
  if isinstance(range_positions[0], list):
    # check if the value is a list of lists; if so, check both lists
    print("is list")
    return is_within_range(position, range_positions[0]) or is_within_range(position, range_positions[1])

  if len(position) > 1:
    print("mulitple positions")
    # if the value is a list of two items, check if the position is within the range
    if any([x in range(range_positions[0], range_positions[1]) for x in position]):
      return True
    if any([x in range(position[0], position[1]) for x in range_positions]):
      return True
 
  # the position is a single item
  elif range_positions[0] <= position[0] <= range_positions[1]:
    print('single position')
    return True
  
  print('not in positions')
  return False



"""
A dictionary to turn TBProfiler WHO annotations into
their corresponding CDPH Looker or MDL interpretations
"""
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
A dictionary that matches CDPH LIMS antimicrobial codes to the 
corresponding drug name; used to create the LIMS report
"""
global ANTIMICROBIAL_CODE_TO_DRUG_NAME
ANTIMICROBIAL_CODE_TO_DRUG_NAME = {
  "M_DST_B01_INH": "isoniazid", 
  "M_DST_C01_ETO": "ethionamide",
  "M_DST_D01_RIF": "rifampin", 
  "M_DST_E01_PZA": "pyrazinamide",
  "M_DST_F01_EMB": "ethambutol",
  "M_DST_G01_AMK": "amikacin", 
  "M_DST_H01_KAN": "kanamycin",
  "M_DST_I01_CAP": "capreomycin", 
  "M_DST_J01_MFX": "moxifloxacin",
  "M_DST_K01_LFX": "levofloxacin", 
  "M_DST_L01_BDQ": "bedaquiline",
  "M_DST_M01_CFZ": "clofazimine", 
  "M_DST_N01_LZD": "linezolid" 
}

"""
A dictionary that matches the LIMS antimicrobial code to the relevant
antimicrobial drug name for cycloserine; activated by the
--add_cs_lims flag
"""
global ANTIMICROBIAL_CODE_TO_DRUG_NAME_CS
ANTIMICROBIAL_CODE_TO_DRUG_NAME_CS = {
  "M_DST_P01_CS": "cycloserine"
}

"""
The dictionary that matches CDPH LIMS antimicrobial codes to the
corresponding genes and their corresponding CDPH LIMS codes;
used to create the LIMS report
"""
global ANTIMICROBIAL_CODE_TO_GENES
ANTIMICROBIAL_CODE_TO_GENES = {}
  
"""
A dictionary that matches CDPH LIMS antimicrobial codes to the
corresponding genes and their corresponding CDPH LIMS codes for tNGS data;
used to create the LIMS report
"""
global ANTIMICROBIAL_CODE_TO_GENES_tNGS
ANTIMICROBIAL_CODE_TO_GENES_tNGS = {
  "M_DST_B01_INH": {
    "katG": "M_DST_B02_katG", 
    "fabG1": "M_DST_B03_fabG1",
    "inhA": "M_DST_B04_inhA"
  },
  "M_DST_C01_ETO": {
    "ethA": "M_DST_C02_ethA", 
    "fabG1": "M_DST_C03_fabG1",
    "inhA": "M_DST_C04_inhA"
  },
  "M_DST_D01_RIF": {
    "rpoB": "M_DST_D02_rpoB"
  },
  "M_DST_E01_PZA": {
    "pncA": "M_DST_E02_pncA"
  },
  "M_DST_F01_EMB": {
    "embB": "M_DST_F03_embB"
  },
  "M_DST_G01_AMK": {
    "rrs": "M_DST_G02_rrs", 
    "eis": "M_DST_G03_eis"
  },
  "M_DST_H01_KAN": {
    "rrs": "M_DST_H02_rrs", 
    "eis": "M_DST_H03_eis"
  },
  "M_DST_I01_CAP": {
    "rrs": "M_DST_I02_rrs", 
    "tlyA": "M_DST_I03_tlyA"
  },
  "M_DST_J01_MFX": {
    "gyrA": "M_DST_J02_gyrA", 
    "gyrB": "M_DST_J03_gyrB"
  },
  "M_DST_K01_LFX": {
    "gyrA": "M_DST_K02_gyrA", 
    "gyrB": "M_DST_K03_gyrB"
  },
  "M_DST_L01_BDQ": {
    "Rv0678": "M_DST_L02_Rv0678", 
  },
  "M_DST_M01_CFZ": {
    "Rv0678":"M_DST_M02_Rv0678", 
  },
  "M_DST_N01_LZD": {
    "rrl": "M_DST_N02_rrl", 
    "rplC": "M_DST_N03_rplC"
  }
}

"""
A dictionary that matches CDPH LIMS antimicrobial codes to the
corresponding genes and their corresponding CDPH LIMS codes for WGS data;
used to create the LIMS report
"""
global ANTIMICROBIAL_CODE_TO_GENES_WGS
ANTIMICROBIAL_CODE_TO_GENES_WGS = {
  "M_DST_B01_INH": {
    "katG": "M_DST_B02_katG", 
    "fabG1": "M_DST_B03_fabG1",
    "inhA": "M_DST_B04_inhA"
  },
  "M_DST_C01_ETO": {
    "ethA": "M_DST_C02_ethA", 
    "fabG1": "M_DST_C03_fabG1",
    "inhA": "M_DST_C04_inhA"
  },
  "M_DST_D01_RIF": {
    "rpoB": "M_DST_D02_rpoB"
  },
  "M_DST_E01_PZA": {
    "pncA": "M_DST_E02_pncA"
  },
  "M_DST_F01_EMB": {
    "embA": "M_DST_F02_embA", 
    "embB": "M_DST_F03_embB"
  },
  "M_DST_G01_AMK": {
    "rrs": "M_DST_G02_rrs", 
    "eis": "M_DST_G03_eis"
  },
  "M_DST_H01_KAN": {
    "rrs": "M_DST_H02_rrs", 
    "eis": "M_DST_H03_eis"
  },
  "M_DST_I01_CAP": {
    "rrs": "M_DST_I02_rrs", 
    "tlyA": "M_DST_I03_tlyA"
  },
  "M_DST_J01_MFX": {
    "gyrA": "M_DST_J02_gyrA", 
    "gyrB": "M_DST_J03_gyrB"
  },
  "M_DST_K01_LFX": {
    "gyrA": "M_DST_K02_gyrA", 
    "gyrB": "M_DST_K03_gyrB"
  },
  "M_DST_L01_BDQ": {
    "Rv0678": "M_DST_L02_Rv0678", 
    "atpE": "M_DST_L03_atpE",
    "pepQ": "M_DST_L04_pepQ", 
    "mmpL5": "M_DST_L05_mmpL5",
    "mmpS5": "M_DST_L06_mmpS5"
  },
  "M_DST_M01_CFZ": {
    "Rv0678":"M_DST_M02_Rv0678", 
    "pepQ": "M_DST_M03_pepQ",
    "mmpL5":"M_DST_M04_mmpL5", 
    "mmpS5": "M_DST_M05_mmpS5"
  },
  "M_DST_N01_LZD": {
    "rrl": "M_DST_N02_rrl", 
    "rplC": "M_DST_N03_rplC"
  }
}

"""
Additional optional entries for the LIMS report for cycloserine
Activated by the --add_cs_lims flag
"""
global ANTIMICROBIAL_CODE_TO_GENES_CS
ANTIMICROBIAL_CODE_TO_GENES_CS = {
  "M_DST_P01_CS": {
  "ald": "M_DST_P02_ald",
  "alr": "M_DST_PO3_alr"
  }
}

"""
A dictionary that contains a list of the antimicrobial drug names;
used to create the Looker report
"""
global ANTIMICROBIAL_DRUG_NAME_LIST
ANTIMICROBIAL_DRUG_NAME_LIST = [
  "amikacin", "bedaquiline", "capreomycin", "clofazimine", 
  "ethambutol", "ethionamide", "isoniazid", "kanamycin", 
  "levofloxacin", "linezolid", "moxifloxacin", "pyrazinamide", 
  "rifampin", "streptomycin"
]

"""
A dictionary that matches the antimicrobial drugs to the
genes that may confer resistance to them; 
see also https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv;
used to create the Looker report
"""
global ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME
ANTIMICROBIAL_DRUG_NAME_TO_GENE_NAME = { 
  "amikacin": ["eis", "rrs"], 
  "bedaquiline": ["Rv0678"], 
  "capreomycin": ["rrs", "tlyA"], 
  "clofazimine": ["Rv0678"], 
  "ethambutol": ["embA", "embB", "embC", "embR"], 
  "ethionamide": ["ethA", "ethR", "fabG1", "inhA"], 
  "isoniazid": ["ahpC", "fabG1", "inhA", "kasA", "katG"], 
  "kanamycin": ["eis", "rrs"], 
  "levofloxacin": ["gyrA", "gyrB"], 
  "linezolid": ["rplC", "rrl"],
  "moxifloxacin": ["gyrA", "gyrB"],
  "pyrazinamide": ["panD", "pncA", "rpsA"], 
  "rifampin": ["rpoA", "rpoB", "rpoC"], 
  "streptomycin": ["gid", "rpsL", "rrs"] 
}

"""
A dictionary that will contain the percent coverage for each gene
"""
global COVERAGE_DICTIONARY
COVERAGE_DICTIONARY = {}
 
"""
The coverage threshold (the minimum percent coverage of a gene/locus over the minimum depth)
See also MIN_DEPTH
"""
global COVERAGE_THRESHOLD
COVERAGE_THRESHOLD = 100

"""
This dataframe is made to host the laboratorian report
"""
global DF_LABORATORIAN
DF_LABORATORIAN = pd.DataFrame(columns = [
  "sample_id", "tbprofiler_gene_name", "tbprofiler_locus_tag", "tbprofiler_variant_substitution_type", 
  "tbprofiler_variant_substitution_nt", "tbprofiler_variant_substitution_aa", "confidence", "antimicrobial",
  "looker_interpretation", "mdl_interpretation", "depth", "frequency", "read_support", "rationale", "warning"
])

"""
The minimum frequency for a mutation to be considered for ethA at protein position 237
"""
global ETHA237_FREQUENCY
ETHA237_FREQUENCY = 0.1

"""
A list of genes that correspond to a certain set of expert rules
Rv0678 is equivalent to mmpR5
"""
global GENE_LIST 
GENE_LIST = ["atpE", "mmpL5", "mmpS5", "pepQ", "rplC", "rrl", "Rv0678", "ethA", "gid", "katG", "pncA", "rpoB"]

"""
A list of genes that correspond to rule 1.1
Rv0678 is equivalent to mmpR5
"""
global GENE_LIST_MDL_1_1
GENE_LIST_MDL_1_1 = ["atpE", "mmpL5", "mmpS5", "pepQ", "rplC", "rrl", "Rv0678"]

"""
A list of genes that correspond to rule 2.1
Rv0678 is equivalent to mmpR5
"""
global GENE_LIST_MDL_2_1
GENE_LIST_MDL_2_1 = ["ethA", "gid", "katG", "pncA", "rpoB"]

"""
A dictionary corresponding each gene to the drug they may confer
resistance to, including the genes in the TBDB watchlist.
See also: https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv, and
https://github.com/jodyphelan/tbdb/blob/master/tbdb.watchlist.csv;
used to create the Laboratorian report
"""
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
  "gyrA": ["ciprofloxacin", "fluoroquinolones", "levofloxacin", "moxifloxacin", "ofloxacin"],
  "gyrB": ["ciprofloxacin", "fluoroquinolones", "levofloxacin", "moxifloxacin", "ofloxacin"],
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
  "rrs": ["amikacin", "aminoglycosides", "capreomycin", "kanamycin", "linezolid"],
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
A dictionary that matches each gene to its corresponding locus tag;
see also: https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb.bed;
used to create the Laboratorian report
"""  
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
A dictionary that matches each gene to it's corresponding tier;
tier information provided by CDPH; used to create the Laboratorian report
"""   
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
The list of genes used to generate the CDPH LIMS report;
used to create the LIMS report
"""
global GENES_FOR_LIMS
GENES_FOR_LIMS = []

"""
A list of genes that will be included in the LIMS report
if --add_cs_lims is true (cycloserine)
"""
global GENES_FOR_LIMS_CS
GENES_FOR_LIMS_CS = [
  "ald", "alr"
]

"""
A  list of genes that are to be included in the 
CDPH LIMS report for tNGS data; used to create the LIMS report
"""
global GENES_FOR_LIMS_tNGS
GENES_FOR_LIMS_tNGS = [
  "eis", "embB", "ethA", "fabG1", "gyrA", 
  "gyrB", "inhA", "katG", "pncA", "rplC", 
  "rpoB", "rrl", "rrs", "Rv0678", "tlyA"
]

"""
A list of genes that are to be included in the 
CDPH LIMS report for WGS data; used to create the LIMS report
"""
global GENES_FOR_LIMS_WGS
GENES_FOR_LIMS_WGS = [
  "atpE", "eis", "embA", "embB", "ethA", "fabG1", "gyrA", 
  "gyrB", "inhA", "katG", "mmpL5", "mmpS5", "pepQ", 
  "pncA", "rplC", "rpoB", "rrl", "rrs", "Rv0678", "tlyA"
]

"""
A set of genes that have been reported so far
"""
global GENES_REPORTED
GENES_REPORTED = set()
  
"""
A set that contains genes that have deletions
"""
global GENES_WITH_DELETIONS
GENES_WITH_DELETIONS = set()

"""
A string that contains the lineage of the sample from the "main_lineage" field
"""
global LINEAGE
LINEAGE = ""

"""
A string that contains the lineage of the sample in English
"""
global LINEAGE_ENGLISH
LINEAGE_ENGLISH = ""
            
                
"""
A list that will contain the names of genes that have coverages
below the coverage threshold (see also COVERAGE_THRESHOLD)
"""  
global LOW_DEPTH_OF_COVERAGE_LIST
LOW_DEPTH_OF_COVERAGE_LIST = []

"""
The minimum depth threshold
"""
global MIN_DEPTH
MIN_DEPTH = 10

"""
The minimum frequency for a mutation to pass QC
"""
global MIN_FREQUENCY
MIN_FREQUENCY = 0.1

"""
The minimum read support for a mutation to pass QC
(calculated as DEPTH * FREQUENCY)
"""
global MIN_READ_SUPPORT
MIN_READ_SUPPORT = 10

""" 
A list that will contain the nucleotide mutations that have failed
QC: min coverage < 10x, min freq < 10%, min read support < 10x (or whatever minimums the user specifies)
"""
global MUTATION_FAIL_LIST
MUTATION_FAIL_LIST = []

""" 
This variable holds the operator's name
"""
global OPERATOR
OPERATOR = ""
          
"""
A dictionary that matches certain genes to their promoter regions.
If a mutation is within these promoter regions, it needs special consideration; these are nucleotide positions
DEPRECATED as of v2.0.1
"""
global PROMOTER_REGIONS
PROMOTER_REGIONS = {
  "atpE": [-48, -1], # CDPH range
  "pepQ": [-33, -1], # CDPH range
  "rplC": [-18, -1], # CDPH range
  "Rv0678": [-84, -1], # CDPH range
}

"""
A dictionary that matches the promoter regions from the WHO v2 catalogue.
If a mutation is within these regions, it needs special consideration; these are nucleotide positions
"""
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
  "whiB7": [-404, -1],
}

"""
A dictionary that ranks the resistances on severity
"""
global RESISTANCE_RANKING
RESISTANCE_RANKING = {
  "R": 6,
  "R-Interim": 5,
  "U": 4,
  "S-Interim": 3,
  "S": 2,
  "WT": 1,
  "Insufficient Coverage": 0,
  "NA": 0
}

"""
A list of rpoB mutations that require unique LIMS output wording
"""
global RPOB_MUTATIONS
RPOB_MUTATIONS = [
  "Leu430Pro",
  "Asp435Tyr",
  "His445Asn",
  "His445Cys",
  "His445Leu",
  "His445Ser",
  "Leu452Pro",
  "Ile491Phe"
]

"""
The minimum frequency for a mutation to be considered for rpoB at protein position 449
"""
global RPOB449_FREQUENCY
RPOB449_FREQUENCY = 0.1

"""
The minimum frequency for a mutation to be considered for rrl
"""
global RRL_FREQUENCY
RRL_FREQUENCY = 0.1

"""
The minimum read support for a mutation to be considered for rrl
"""
global RRL_READ_SUPPORT
RRL_READ_SUPPORT = 10
            
"""
The minimum frequency for a mutation to be considered for rrs
"""            
global RRS_FREQUENCY
RRS_FREQUENCY = 0.1

"""
The minimum read support for a mutation to be considered for rrs
"""
global RRS_READ_SUPPORT
RRS_READ_SUPPORT = 10
            
"""
This dictionary converts the rule to the rationale language
"""
global RULE_TO_RATIONALE
RULE_TO_RATIONALE = {
  "rule1.2": "Expert rule 1.2. Novel drug targets",
  "rule2.2.1": "Expert rule 2.2.1. Loss-of-function",
  "rule2.2.2.1": "Expert rule 2.2.2.1. rpoB RRDR",
  "rule2.2.2.2": "Expert rule 2.2.2.2 rpoB non-RRDR",
  "rule3.2.1": "Expert rule 3.2.1. rrs",
  "rule3.2.2": "Expert rule 3.2.3. gyrA QRDR",
  "rule3.2.3": "Expert rule 3.2.3. gyrB QRDR",
  "rule3.2.4": "No WHO annotation or expert rule",
  "rule6": "Mutation in proximal promoter region"
}

"""
This variable holds the sample name
"""
global SAMPLE_NAME
SAMPLE_NAME = ""

"""
This variable holds the sequencing method
"""
global SEQUENCING_METHOD
SEQUENCING_METHOD = ""

"""
This is a dictionary of positions for genes requiring different consideration.
Note: the rpoB, gyrA, and gyrB special positions are in codons, rrl & rrs are nucleotide positions
rpoB and rrl indicate ranges; rrs indicates specific positions
"""
global SPECIAL_POSITIONS
SPECIAL_POSITIONS = {
  "rpoB": [426, 452], # codon; the RRDR range
  "gyrA": [88, 94], # codon; the QRDR range
  "gyrB": [446, 507], # codon; the QRDR range
  "rrl": [[2003, 2367], [2449, 3056]], # nucleotide; range
  "rrs": [1401, 1402, 1484] # nucleotide; specific positions
}

"""
The specific primer regions for the gene(s) that is(are) split into multiple sections [tNGS only]
Left blank to not affect WGS sequencing
"""
global TNGS_REGIONS
TNGS_REGIONS = {}

"""
The specific primer regions for the gene(s) that is(are) split into multiple sections [tNGS only]
"""
global TNGS_REGIONS_ACTIVATED
TNGS_REGIONS_ACTIVATED = {
  "rpoB": {
    "rpoB_1": [760957, 761355],
    "rpoB_2": [760280, 760812]
  }
}