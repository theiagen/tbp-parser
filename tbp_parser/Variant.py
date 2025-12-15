import re
import globals as globals_
from Row import Row 
import copy

class Variant:
    """
    This class represents a single Variant reported by TBProfiler.
    
    Note:
    A variant can have either (1) an annotation field. This field could
    (a) have 1 or more annotations, or (b) have none. A variant could also 
    have (2) no annotation field at all.

    Within the annotation field, an annotation could have more than one drug
    listed. In addition, it is possible that the same drug can show up twice
    within the same annotation for a single variant.

    This class has three functions:
        - extract_annotations: splits up the annotation field into individual
            annotations and creates a Row object for each annotation.
        - extract_consequences: iterates through the consequences section of
            the variant to add any additional genes that this mutation could be
            classified under.
        - apply_expert_rules: applies the expert rules from the interpretation
            logic document to the variant.
    """
    def __init__(self, logger, sample_name, gene_to_locus_tag, gene_to_tier, gene_to_antimicrobial_drug_name, promoter_regions, variant=None):
        """Initializes the Variant class

        Args:
            logger (logging): generates logs
            sample_name (str): the sample name
            gene_to_locus_tag (dict): a dictionary mapping genes to their locus tags
            gene_to_tier (dict): a dictionary mapping genes to their tiers
            gene_to_antimicrobial_drug_name (dict): a dictionary mapping genes to their associated antimicrobial drugs
            variant (dict, optional): a dictionary of attributes about the variant. Defaults to None. The potential keys in the dictionary are:

            - chrom (str): the name of the chromosome in the BAM file
            - pos (int): the position of the variant (in reference to the entire genome)
            - ref (str): the reference nucleotide(s)
            - alt (str): the alternate nucleotide(s)
            - depth (int): the depth of coverage at this position
            - freq (float): the frequency of the variant at this position
            - sv (bool): whether the variant is a structural variant
            - filter (str): filter status of the variant
            - forward_reads (int): number of forward reads supporting the variant
            - reverse_reads (int): number of reverse reads supporting the variant
            - sv_len (int): length of the structural variant; null if not a structural variant
            - gene_id (str): the gene ID where the variant is located
            - gene_name (str): the gene name where the variant is located
            - feature_id (str): the feature ID where the variant is located
            - type (str): the type of variant (e.g., synonymous_variant)
            - change (str): the protein change caused by the variant, or the nucleotide change if non-coding
            - nucleotide_change (str): the nucleotide change caused by the variant
            - protein_change (str): the protein change caused by the variant
            - annotation (list): a list of annotations associated with the variant
            - consequences (list): a list of alternate consequences associated with the variant
            - drugs (list): a list of drugs associated with the gene and the resistance call and/or comment for each drug
            - locus_tag (str): the locus tag of the gene where the variant is located
            - gene_associated_drugs (list): a list of drugs associated with the gene where the variant is located
        """        
        self.logger = logger

        self.annotation_dictionary = {}
        """A dictionary containing the various annotations for this variant"""
        
        self.sample_name = sample_name
        """The sample name associated with this variant."""
        
        self.GENE_TO_LOCUS_TAG = gene_to_locus_tag
        """A dictionary mapping genes to their locus tags."""

        self.GENE_TO_TIER = gene_to_tier
        """A dictionary mapping genes to their tiers."""
        
        self.GENE_TO_ANTIMICROBIAL_DRUG_NAME = gene_to_antimicrobial_drug_name
        """A dictionary mapping genes to their associated antimicrobial drugs."""

        self.PROMOTER_REGIONS = promoter_regions
        """A dictionary mapping genes to their promoter regions."""

        if variant is not None:
            for key, value in variant.items():
                # see docstring for possible keys
                setattr(self, key, value)

    def extract_annotations(self) -> None:
        """
        This function takes the annotation field in a variant and splits it into
        its individual parts, creating an Row object for the annotation field. If the
        annotation field is empty or does not exist, the function creates a row
        based off of the gene_associated_drugs field.
        """
        mock_annotation = {
            "confidence": "No WHO annotation",
            "drug": "",
            "source": "",
            "comment": ""
        }
        
        if hasattr(self, "annotation") and len(self.annotation) > 0:
            # if possibility 1a (variant has an annotation field with content (len > 0))
            # create a list of the drugs associated with the gene to check if all drugs are reported
            gene_associated_drug_list = []
            if hasattr(self, "gene_associated_drugs"):
                gene_associated_drug_list = self.gene_associated_drugs
        
            # iterate through each entry in the annotation
            for annotation_entry in self.annotation:
                drug = annotation_entry["drug"]
                
                new_row = Row(self.logger, self, annotation_entry)
                
                # if this is the first time a drug has been seen, add it to the annotation dictionary
                if drug not in self.annotation_dictionary.keys():
                    self.annotation_dictionary[drug] = new_row

                    # if this drug is in the gene associated drug list, remove it because it has been accounted for
                    if drug in gene_associated_drug_list:
                        gene_associated_drug_list.remove(drug)

                # if it has been seen before, save the row with the more severe WHO confidence (higher value)
                elif new_row.rank_annotation() > self.annotation_dictionary[drug].rank_annotation():
                    self.annotation_dictionary[drug] = new_row
                # also preferentially keep WHO confidence rows over non-WHO confidence rows if the ranks are the same    
                elif new_row.rank_annotation() == self.annotation_dictionary[drug].rank_annotation():
                    if "WHO" in new_row.source and "WHO" not in self.annotation_dictionary[drug].source:
                        self.annotation_dictionary[drug] = new_row
         
            # add any missing drugs from the annotation list that are found on the gene associated drug list to the annotation dictionary 
            for drug in gene_associated_drug_list:  
                # see also rule 4.3.3
                # copy a row in the annotation dictionary, but update the confidence and drug name appropriately. no qc has been performed so this should be fine
                mock_annotation["drug"] = drug
                
                self.annotation_dictionary[drug] = Row(self.logger, self, mock_annotation)

            if self.gene_name in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME.keys():
                for drug in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME[self.gene_name]:
                    if drug not in self.annotation_dictionary.keys():
                        # copy a row in the annotation dictionary, but update the confidence, source and drug name appropriately. no qc has been performed so this should be fine
                        mock_annotation["source"] = "Mutation effect for given drug is not in TBDB"
                        mock_annotation["drug"] = drug
                                                
                        self.annotation_dictionary[drug] = Row(self.logger, self, mock_annotation)

        else:
            # possibilities 1b and 2: the annotation field has no content or the field does not exist
            for drug in self.gene_associated_drugs:
                mock_annotation["drug"] = drug
                
                self.annotation_dictionary[drug] = Row(self.logger, self, mock_annotation)

            if self.gene_name in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME.keys():
                for drug in self.GENE_TO_ANTIMICROBIAL_DRUG_NAME[self.gene_name]:
                    if drug not in self.annotation_dictionary.keys():
                        mock_annotation["drug"] = drug
                        
                        self.annotation_dictionary[drug] = Row(self.logger, self, mock_annotation)

    def extract_consequences(self) -> set:
        """This section iterates through the "consequences" section of a variant in 
        order to add any additional genes that this mutation could be classified under.
        For example, mutations in mmpR5 could also be classified under mmpS5 or mmpL5.
        
        Returns:
            set: a set of gene names that were reported from the consequences section
        """        
        reported_genes = set()
    
        self.logger.debug("VAR:extract_consequences:The gene is {}, now checking for consequences".format(self.gene_name))
        if hasattr(self, "consequences") and len(self.consequences) > 0:
            # if there is a consequences section that has content, iterate
            for entry in self.consequences:
                if entry.get("gene_name") != self.gene_name:
                    # make a copy of the variant object 
                    copied_variant = copy.deepcopy(self) 
                    
                    # remove drugs, locus_tag, gene_associated_drugs, and the annotation field from the copied variant
                    attrs_to_remove = ["change", "annotation", "consequences", "drugs", "locus_tag", "gene_associated_drugs"]
                    for attr in attrs_to_remove:
                        if hasattr(copied_variant, attr):
                            delattr(copied_variant, attr)

                    for key, value in entry.items():
                        # reset the attributes of the copied variant to those of the consequence entry
                        setattr(copied_variant, key, value)
                    
                    copied_variant.extract_annotations()
                    
                    reported_genes.add(copied_variant.gene_name)
                    
        return reported_genes 

    def apply_expert_rules(self) -> tuple[str, str]:
        """Applies the expert rules from the interpretation logic document regarding 
        the interpretation of potential resistance mutations in the absence of WHO 
        annotations.

        Returns:
            tuple[str, str]: The interpretation (str) and the rule applied (str)
        """        
        GENE_LIST = ["atpE", "mmpL5", "mmpS5", "pepQ", "rplC", "rrl", "Rv0678", "ethA", "gid", "katG", "pncA", "rpoB"]
        """A list of genes that correspond to a certain set of expert rules."""

        SPECIAL_POSITIONS = {
            "rpoB": [426, 452], # codon; the RRDR range
            "gyrA": [88, 94], # codon; the QRDR range
            "gyrB": [446, 507], # codon; the QRDR range
            "rrl": [[2003, 2367], [2449, 3056]], # nucleotide; range
            "rrs": [1401, 1402, 1484] # nucleotide; specific positions
        }
        """This is a dictionary of special positions for genes requiring different consideration.
        Note: the rpoB, gyrA, and gyrB special positions are in codons, 
        rrl & rrs are nucleotide positions
        all are ranges except rrs indicates specific positions
        """
        
        position_nt = globals_.get_position(self.nucleotide_change)
        position_aa = globals_.get_position(self.protein_change)

        self.logger.debug("VAR:apply_expert_rules:The nucleotide position is {} and the amino acid position is {}".format(position_nt, position_aa))
        
        # FYI: There is a lot of duplicated code here but it is kept this way to ensure
        # that the rule is appropriately assigned based on the gene and scenario
        if self.gene_name in ["Rv0678", "atpE", "pepQ", "rplC", "mmpL5", "mmpS5"]:        
            self.logger.debug("VAR:apply_expert_rules:The gene is {}, now checking if the position requires special consideration under rule 1.2".format(self.gene_name))

            # check if position within promoter regions
            if self.gene_name in self.PROMOTER_REGIONS.keys() and globals_.is_within_range(position_nt, self.PROMOTER_REGIONS[self.gene_name]):
                self.logger.debug("VAR:apply_expert_rules:The position is within the promoter region; interpretation is 'U'")
                return "U", "whov2"

            # otherwise, check if it is an upstream gene variant
            if "upstream_gene_variant" in self.type: 
                self.logger.debug("VAR:apply_expert_rules:The position is an upstream gene variant and NOT in the proximal promoter regions; interpretation is 'S'")
                return "S", "rule1.2"

            # otherwise, check if it is not in the ORF
            if ((not any(non_ORF in self.nucleotide_change for non_ORF in ["+", "-", "*"]) or self.nucleotide_change.endswith("*")) 
                or (not any(non_ORF in self.protein_change for non_ORF in ["+", "-", "*"]) or self.protein_change.endswith("*"))):
                self.logger.debug("VAR:apply_expert_rules:The position is not in the ORF; interpretation is 'S' if it is a synonymous variant or 'U' if it is not")
                # if a position includes either +, *, or - it's not in the ORF, unless 
                #  the * is at the end which means its a premature stop codon
                if self.type != "synonymous_variant":
                    return "U", "rule1.2"
                else:
                    return "S", "rule1.2"

        elif self.gene_name == "rrl":
            self.logger.debug("VAR:apply_expert_rules:The gene is rrl, now checking if the position requires special consideration")

            if globals_.is_within_range(position_nt, SPECIAL_POSITIONS[self.gene_name]):
                return "U", "rule1.2"

            # check if position within promoter regions
            if self.gene_name in self.PROMOTER_REGIONS.keys() and globals_.is_within_range(position_nt, self.PROMOTER_REGIONS[self.gene_name]):
                self.logger.debug("VAR:apply_expert_rules:The position is within the proximal promoter region; interpretation is 'U'")
                return "U", "whov2"

            else:
                self.logger.debug("VAR:apply_expert_rules:The position is not within the special positions or the proximal promoter region; interpretation is 'S'")   
                return "S", "rule1.2"

        elif self.gene_name in ["katG", "pncA", "ethA", "gid"]: 
            self.logger.debug("VAR:apply_expert_rules:The gene is {}, now checking if the mutation type requires special consideration under rule 2.2".format(self.gene_name))

            if ((any(indel_or_stop in self.nucleotide_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or self.nucleotide_change.endswith("*")) 
                or (any(indel_or_stop in self.protein_change for indel_or_stop in ["del", "ins", "fs", "delins", "_"]) or self.protein_change.endswith("*"))):
                if any([int(position) > -30 for position in position_nt]): 
                    self.logger.debug("VAR:apply_expert_rules:The mutation type is an indel, stop, or frameshift codon and within 30 nt of the start codon; interpretation is 'R'")
                    return "R", "rule2.2.1"
                else:
                    self.logger.debug("VAR:apply_expert_rules:The mutation type is an indel, stop, or frameshift codon but is past 30 nt of the start codon; interpretation depends on variant type")

            if (self.type == "synonymous_variant"):
                self.logger.debug("VAR:apply_expert_rules:The mutation type is a synonymous variant; interpretation is 'S'")
                return "S", "rule2.2.1"

            # check if position within promoter regions
            if self.gene_name in self.PROMOTER_REGIONS.keys() and globals_.is_within_range(position_nt, self.PROMOTER_REGIONS[self.gene_name]):
                self.logger.debug("VAR:apply_expert_rules:The position is within the proximal promoter region; interpretation is 'U'")
                return "U", "whov2"

            elif ("upstream_gene_variant" in self.type):
                self.logger.debug("VAR:apply_expert_rules:The mutation type is a NONsynonymous variant, not in the proximal promoter region, and is an upstream gene variant; interpretation is 'S'")     
                return "S", "rule2.2.1"
            else:
                self.logger.debug("VAR:apply_expert_rules:The mutation type IS a NONsynonymous variant, not in the proximal promoter region, and is NOT an upstream gene variant; interpretation is 'U'")
                return "U", "rule2.2.1"

        # rules 2.2.2.1 and 3.2.2 & 3.2.3
        elif self.gene_name in ["gyrA", "gyrB", "rpoB"]: 
            self.logger.debug("VAR:apply_expert_rules:The gene is {}, now checking if the position requires special consideration".format(self.gene_name))

            if globals_.is_within_range(position_aa, SPECIAL_POSITIONS[self.gene_name]):
                self.logger.debug("VAR:apply_expert_rules:The position is within the special positions; interpretation is 'R' if rpoB (or 'U' if not) and nonsynonymous, else 'S'")

                if self.gene_name == "rpoB":
                    if self.type == "synonymous_variant":
                        return "S", "rule2.2.2.1" 
                    else:
                        return "R", "rule2.2.2.1"

                if self.gene_name == "gyrA":
                    if self.type != "synonymous_variant":
                        return "U", "rule3.2.2"

                if self.gene_name == "gyrB":
                    if self.type != "synonymous_variant":
                        return "U", "rule3.2.3"

            elif (self.type == "synonymous_variant"):
                self.logger.debug("VAR:apply_expert_rules:The position is not within the special positions and is synonymous; interpretation is 'S'")
                return ("S", "rule2.2.2.2") if self.gene_name == "rpoB" else ("S", "rule3.2.4")

            if self.gene_name in self.PROMOTER_REGIONS.keys() and globals_.is_within_range(position_nt, self.PROMOTER_REGIONS[self.gene_name]):
                self.logger.debug("VAR:apply_expert_rules:The position is within the proximal promoter region; interpretation is 'U'")
                return "U", "whov2"

            elif ("upstream_gene_variant" in self.type):
                self.logger.debug("VAR:apply_expert_rules:The position is not within the special positions, not in the proximal promoter region, is nonsynomyous but is an upstream gene variant; interpretation is 'S'")

                if self.gene_name == "rpoB":
                    return "S", "rule2.2.2.2"
                else:
                    return "S", "rule3.2.4"
            else:
                self.logger.debug("VAR:apply_expert_rules:The position is not within the special positions, not in the proximal promoter region, is nonsynonymous and is NOT an upstream gene variant; interpretation is 'U'")
                return ("U", "rule2.2.2.2") if self.gene_name == "rpoB" else ("U", "rule3.2.4")

        elif self.gene_name not in GENE_LIST:
            self.logger.debug("VAR:apply_expert_rules:The gene is not in the gene list that requires an expert rule.")

            if self.gene_name == "rrs":
                self.logger.debug("VAR:apply_expert_rules:The gene is rrs, now checking if the position requires special consideration")

                if any(map(lambda position: position in position_nt, SPECIAL_POSITIONS[self.gene_name])):
                    self.logger.debug("VAR:apply_expert_rules:The position is within the special positions; interpretation is 'U'")
                    return "U", "rule3.2.1"

                elif self.gene_name in self.PROMOTER_REGIONS.keys() and globals_.is_within_range(position_nt, self.PROMOTER_REGIONS[self.gene_name]):
                    self.logger.debug("VAR:apply_expert_rules:The position is within the proximal promoter region; interpretation is 'U'")
                    return "U", "whov2"

                else:
                    self.logger.debug("VAR:apply_expert_rules:The position is not within the special positions, and not in the proximal promoter region,; interpretation is 'S'")
                    return "S", "rule3.2.1"

            # rule 3.2.4: all remaining scenarios not covered above
            elif (self.type == "synonymous_variant"):
                self.logger.debug("VAR:apply_expert_rules:The mutation is synonymous; interpretation is 'S'")
                return "S", "rule3.2.4"
            elif self.gene_name in self.PROMOTER_REGIONS.keys() and globals_.is_within_range(position_nt, self.PROMOTER_REGIONS[self.gene_name]):
                self.logger.debug("VAR:apply_expert_rules:The position is within the proximal promoter region; interpretation is 'U'")
                return "U", "whov2"
            elif ("upstream_gene_variant" in self.type):
                self.logger.debug("VAR:apply_expert_rules:The mutation is a nonsynonymous variant and not in the proximal promoter region, but is an upstream gene variant; interpretation is 'S'")
                return "S", "rule3.2.4"
            else:
                self.logger.debug("VAR:apply_expert_rules:The mutation is a nonsynonymous variant, not in the proximal promoter region, and is NOT an upstream gene variant; interpretation is 'U'")
                return "U", "rule3.2.4"

        return "", ""
    