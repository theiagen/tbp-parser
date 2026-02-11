import json
import logging
import pandas as pd
from datetime import datetime
from pathlib import Path

from variant import Variant
from utils.config import Configuration
from utils.helper import Helper

logger = logging.getLogger(__name__)

# LIMS uses a simplified ranking (no interim designations)
RESISTANCE_RANKING = {
    "R": 4,
    "U": 3,
    "S": 2,
    "WT": 1,
    "Insufficient Coverage": 0,
    "NA": -1,
}

# Codon/nucleotide positions for genes requiring special LIMS handling
SPECIAL_POSITIONS = {
    "rpoB": [426, 452],
    "gyrA": [88, 94],
    "gyrB": [446, 507],
    "rrl": [[2003, 2367], [2449, 3056]],
    "rrs": [1401, 1402, 1484],
}

# rpoB mutations indicating low-level resistance
RPOB_MUTATIONS = [
    "Leu430Pro", "Asp435Tyr", "His445Asn", "His445Cys",
    "His445Leu", "His445Ser", "Leu452Pro", "Ile491Phe",
]


def _convert_annotation(annotation: str, drug: str) -> str:
    """Convert a resistance annotation into LIMS text."""
    if annotation == "R":
        return f"Mutation(s) associated with resistance to {drug} detected"
    elif annotation in ("R-Interim", "U"):
        return f"The detected mutation(s) have uncertain significance. Resistance to {drug} cannot be ruled out"
    elif annotation == "Insufficient Coverage":
        return "Pending Retest"
    return f"No mutations associated with resistance to {drug} detected"


def get_lineage_id(
    config: Configuration,
    lims_format: dict,
    low_depth_genes: list[str],
    genes_with_valid_deletions: list[str],
) -> tuple[str, str]:
    """Determine the lineage and English lineage string from TBProfiler JSON.

    Args:
        config: Configuration object
        lims_format: LIMS report format dict loaded from YAML
        low_depth_genes: Genes that failed locus coverage QC
        genes_with_valid_deletions: Genes with QC-passing deletions

    Returns:
        tuple[str, str]: (raw lineage from TBProfiler, English lineage string)
    """
    with open(config.input_json) as f:
        input_json = json.load(f)

    detected_lineage = input_json.get("main_lineage", "")
    detected_sublineage = input_json.get("sub_lineage", "")

    # collect all unique gene names from LIMS format
    lims_genes = list({
        gene_name
        for code_to_genes in lims_format.values()
        for gene_to_code in code_to_genes.values()
        for gene_name in gene_to_code.keys()
    })

    # count passing genes
    passing_gene_count = sum(
        1 for gene in lims_genes
        if gene not in low_depth_genes and gene not in genes_with_valid_deletions
    )

    try:
        pct_above = passing_gene_count / len(lims_genes)
    except ZeroDivisionError:
        pct_above = 0

    lineage = set()

    if pct_above >= config.MIN_PERCENT_LOCI_COVERED:
        if config.TNGS:
            # tNGS lineage determination based on pncA His57Asp
            # TODO: port full tNGS lineage logic (pncA His57Asp check with QC)
            lineage.add("DNA of Mycobacterium tuberculosis complex detected (not M. bovis)")
        else:
            sublineages = detected_sublineage.split(";")
            if "lineage" in detected_lineage:
                lineage.add("DNA of Mycobacterium tuberculosis species detected")

            for sublineage in sublineages:
                if "BCG" in detected_lineage or "BCG" in sublineage:
                    lineage.add("DNA of Mycobacterium bovis BCG detected")
                elif ("La1" in detected_lineage or "La1" in sublineage) or \
                     ("bovis" in detected_lineage or "bovis" in sublineage):
                    lineage.add("DNA of Mycobacterium bovis (not BCG) detected")

            if detected_lineage == "" or detected_lineage == "NA" or len(lineage) == 0:
                logger.debug("No lineage detected by TBProfiler; assuming M.tb")
                lineage.add("DNA of Mycobacterium tuberculosis complex detected")
    else:
        lineage.add("DNA of Mycobacterium tuberculosis complex NOT detected")

    lineage_english = "; ".join(sorted(lineage))
    return detected_lineage, lineage_english


def write_lims_report(
    config: Configuration,
    variants: list[Variant],
    low_depth_genes: list[str],
    genes_with_valid_deletions: list[str],
) -> tuple[Path, str, str]:
    """Write the LIMS report from processed Variant objects.

    Args:
        config: Configuration object
        variants: All processed Variant objects (including WT)
        low_depth_genes: Genes that failed locus coverage QC
        genes_with_valid_deletions: Genes with QC-passing deletions

    Returns:
        tuple[Path, str, str]: (path to CSV, raw lineage, English lineage)
    """
    lims_format = config.load_lims_report_format()

    raw_lineage, lineage_english = get_lineage_id(
        config, lims_format, low_depth_genes, genes_with_valid_deletions
    )

    # Build set of variants that failed positional QC (to exclude from LIMS)
    positional_qc_fails = {}
    for v in variants:
        if v.fails_qc:
            if v.gene_name not in positional_qc_fails:
                positional_qc_fails[v.gene_name] = set()
            positional_qc_fails[v.gene_name].add(v.nucleotide_change)

    # Filter to QC-passing variants
    qc_pass_variants = [
        v for v in variants
        if v.nucleotide_change not in positional_qc_fails.get(v.gene_name, set())
    ]

    lims_report = {
        # "Sample_Name": variants[0].sample_id if variants else "",
        # "Lineage_ID": lineage_english,
        "MDL sample accession numbers": variants[0].sample_id if variants else "", # Temp change
        "M_DST_A01_ID": lineage_english, # Temp change
    }

    for drug, drug_gene_dict in lims_format.items():
        for antimicrobial_code, gene_codes in drug_gene_dict.items():
            drug_variants = [
                v for v in qc_pass_variants
                if v.drug == drug and v.gene_name in gene_codes
            ]

            # Determine max MDL resistance
            mdl_interpretations = [v.mdl_interpretation for v in drug_variants if v.mdl_interpretation]
            try:
                max_mdl = max(mdl_interpretations, key=lambda x: RESISTANCE_RANKING.get(x, -1))
            except ValueError:
                max_mdl = "Insufficient Coverage"

            logger.debug(f"LIMS: max MDL resistance for {drug} is {max_mdl}")
            lims_report[antimicrobial_code] = _convert_annotation(max_mdl, drug)

            for gene, gene_code in gene_codes.items():
                gene_variants = [v for v in drug_variants if v.gene_name == gene]
                mutation_list = []

                if max_mdl in ["WT", "Insufficient Coverage", "NA"]:
                    if gene in low_depth_genes and gene not in genes_with_valid_deletions:
                        lims_report[gene_code] = "No sequence"
                        lims_report[antimicrobial_code] = "Pending Retest"
                    else:
                        lims_report[gene_code] = "No mutations detected"

                elif max_mdl == "S":
                    if gene == "rpoB" and drug == "rifampicin":
                        lims_report[antimicrobial_code] = "Predicted susceptibility to rifampicin"
                        synonymous_variants = [v for v in gene_variants if v.type == "synonymous_variant"]
                        for v in synonymous_variants:
                            position_aa = Helper.get_position(v.protein_change)
                            if Helper.is_mutation_within_range(position_aa, SPECIAL_POSITIONS["rpoB"]):
                                lims_report[antimicrobial_code] = (
                                    "Predicted susceptibility to rifampicin. "
                                    "The detected synonymous mutation(s) do not confer resistance"
                                )
                                mutation_list.append(f"{v.protein_change} [synonymous]")
                        lims_report[gene_code] = "; ".join(mutation_list) if mutation_list else "No high confidence mutations detected"
                    else:
                        gene_mdls = [v.mdl_interpretation for v in gene_variants if v.mdl_interpretation]
                        if gene_mdls and max(gene_mdls, key=lambda x: RESISTANCE_RANKING.get(x, -1)) == "S":
                            lims_report[gene_code] = "No high confidence mutations detected"
                        else:
                            lims_report[gene_code] = "No mutations detected"

                elif max_mdl in ["R", "U"]:
                    for v in gene_variants:
                        is_rrdr_syn = (
                            gene == "rpoB"
                            and Helper.is_mutation_within_range(
                                Helper.get_position(v.protein_change),
                                SPECIAL_POSITIONS.get(gene, [])
                            )
                            and v.type == "synonymous_variant"
                            and v.mdl_interpretation == "S"
                        )
                        if v.mdl_interpretation in ["R", "U"] or is_rrdr_syn:
                            if v.type == "synonymous_variant":
                                substitution = f"{v.protein_change} [synonymous]"
                            elif v.protein_change and v.protein_change not in ("", "p.0?", "NA"):
                                substitution = v.protein_change
                            else:
                                substitution = v.nucleotide_change
                            mutation_list.append(substitution)

                    # rpoB special handling for low-level resistance
                    if gene == "rpoB" and max_mdl == "R":
                        rpob_specific_count = 0
                        for mut in mutation_list:
                            if any(rpob_mut in mut for rpob_mut in RPOB_MUTATIONS):
                                rpob_specific_count += 1
                            else:
                                rpob_specific_count = 0
                                break

                        if rpob_specific_count > 0:
                            lims_report[antimicrobial_code] = (
                                "Predicted low-level resistance to rifampicin. "
                                "May test susceptible by phenotypic methods."
                            )
                        else:
                            lims_report[antimicrobial_code] = "Predicted resistance to rifampicin"

                # Write gene_code if mutations found
                if mutation_list:
                    lims_report[gene_code] = "; ".join(mutation_list)
                elif gene_code not in lims_report or lims_report[gene_code] not in ["No sequence", "No mutations detected"]:
                    gene_mdls = [v.mdl_interpretation for v in gene_variants if v.mdl_interpretation]
                    if gene_mdls and max(gene_mdls, key=lambda x: RESISTANCE_RANKING.get(x, -1)) == "S":
                        lims_report[gene_code] = "No high confidence mutations detected"
                    else:
                        lims_report[gene_code] = "No mutations detected"

    # Add metadata
    lims_report["Analysis_Date"] = datetime.now().strftime("%Y-%m-%d %H:%M")
    lims_report["Operator"] = config.OPERATOR

    lims_report["M_DST_O01_Lineage"] = raw_lineage # Temp change
    # lims_report["Lineage"] = raw_lineage

    # Write CSV and transposed CSV
    df = pd.DataFrame([lims_report])
    output_path = Path(f"{config.OUTPUT_PREFIX}.lims_report.csv")
    df.to_csv(output_path, index=False)

    transposed_path = Path(f"{config.OUTPUT_PREFIX}.lims_report.transposed.csv")
    df.T.to_csv(transposed_path, header=False)

    logger.info(f"LIMS report written to {output_path}")
    return output_path, raw_lineage, lineage_english
