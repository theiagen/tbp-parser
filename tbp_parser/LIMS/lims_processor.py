import logging
from typing import Dict, Tuple
from tbp_parser.Coverage.coverage_data import LocusCoverage
from tbp_parser.LIMS.lims_record import LIMSRecord, LIMSGeneCode
from tbp_parser.Utilities.config import Configuration
from tbp_parser.Utilities.gene_database import GeneDatabase
from tbp_parser.Utilities.helper import Helper
from tbp_parser.Variant.variant import Variant

logger = logging.getLogger(__name__)

# LIMS uses a simplified ranking (no interim designations)

class LIMSProcessor:
    """Processor for all decision logic for LIMS report generation (Section 5 rules).
    LIMS report fields based on TBProfiler results and QC data."""

    def __init__(self):
        self.config = Configuration.get_instance()

    RESISTANCE_RANKING = {
        "R": 4,
        "Insufficient Coverage": 3,
        "U": 2,
        "S": 1,
        "WT": 0,
        "NA": -1,
    }

    def process(
        self,
        lims_records: list[LIMSRecord],
        variants: list[Variant],
        locus_coverage_map: Dict[str, LocusCoverage],
        detected_lineage: str,
        detected_sublineage: str,
    ) -> Tuple[list[LIMSRecord], str]:
        """
        Orchestrator function to process LIMS records and determine lineage for LIMS report.
        Args:
            lims_records: List of LIMSRecord objects with populated drug/gene codes but empty interpretations
            variants: List of all Variant objects (including WT)
            locus_coverage_map: Mapping of locus_tag names to their coverage percentages
            detected_lineage: Raw lineage string from TBProfiler JSON
            detected_sublineage: Raw sublineage string from TBProfiler JSON
        Returns:
            - List of LIMSRecord objects with all interpretations populated,
            - Detected lineage for LIMS report
        """
        lims_records = self.process_lims_records(
            lims_records=lims_records,
            locus_coverage_map=locus_coverage_map,
            variants=variants,
        )

        lims_lineage = self.process_lims_mtbc_id(
            lims_records=lims_records,
            variants=variants,
            locus_coverage_map=locus_coverage_map,
            detected_lineage=detected_lineage,
            detected_sublineage=detected_sublineage,
        )
        return lims_records, lims_lineage

    def process_lims_records(
        self,
        lims_records: list[LIMSRecord],
        locus_coverage_map: Dict[str, LocusCoverage],
        variants: list[Variant],
    ) -> list[LIMSRecord]:
        """Process LIMSRecord objects to prepare for LIMS report generation.

        Args:
            lims_records: List of LIMSRecord objects with populated drug/gene codes but empty interpretations
            locus_coverage_map: Mapping of locus_tag names to their coverage percentages
            variants: List of all Variant objects

        Returns:
            List of LIMSRecord objects with all interpretations populated
        """

        for lims_record in lims_records:
            for gene, gene_code in lims_record.gene_codes.items():
                # get candidate/qc_filtered variants for this drug/gene to determine max MDL interpretation
                candidate_variants = [
                    v for v in variants
                    if v.drug == lims_record.drug
                    if v.gene_id == GeneDatabase.get_locus_tag(gene)
                ]

                qc_filtered_variants = []
                for v in candidate_variants:
                    locus_coverage = locus_coverage_map.get(v.gene_id, None)
                    if not locus_coverage:
                        logger.debug(f"No locus coverage found for {v.gene_name}|{v.gene_id} at position {v.pos}")
                        continue

                    # filter out variants that fail positional QC, but keep them if they also fail locus QC
                    if (
                        not v.fails_positional_qc or
                        v.fails_positional_qc and v.fails_locus_qc
                    ):
                        qc_filtered_variants.append(v)

                logger.debug(f"Generating LIMS result - GENE_TARGET - [{lims_record.drug}|{gene}]; {len(qc_filtered_variants)}/{len(candidate_variants)} candidate variants after QC filtering")

                # get max mdl based on filtered variants
                self.set_gene_code_max_mdl(gene_code, qc_filtered_variants)

                # get gene_target_value based on max_mdl interpretation
                self.resolve_gene_target(gene_code)

            self.resolve_drug_target(lims_record)
        return lims_records


    def resolve_drug_target(self, lims_record: LIMSRecord) -> None:
        """Resolve drug target value for a given LIMSRecord based on its gene codes.

        Args:
            lims_record: A LIMSRecord object with populated gene codes and their interpretations

        Returns:
            None
        """
        # rpoB mutations indicating low-level resistance
        RPOB_MUTATIONS = [
            "p.Leu430Pro", "p.Asp435Tyr", "p.His445Asn", "p.His445Cys",
            "p.His445Leu", "p.His445Ser", "p.Leu452Pro", "p.Ile491Phe",
        ]

        all_gene_codes = [
            gene_code for gene_code in lims_record.gene_codes.values()
            if gene_code.max_mdl_interpretation
        ]

        max_gene_code = max(
            all_gene_codes,
            key=lambda gc: self.RESISTANCE_RANKING.get(gc.max_mdl_interpretation or "NA", -1),
        )

        is_rpob_rifampicin = (
            lims_record.drug == "rifampicin" and
            "rpoB" in lims_record.gene_codes and
            max_gene_code == lims_record.gene_codes["rpoB"]
        )
        logger.debug(f"Generating LIMS result - DRUG_TARGET - [{lims_record.drug}|{lims_record.drug_code}] based on results from max gene code: `{max_gene_code.gene_code}`")

        # rifampicin specific drug target logic
        if is_rpob_rifampicin:
            if max_gene_code.max_mdl_interpretation in ["R"]:
                if all([v.protein_change in RPOB_MUTATIONS for v in max_gene_code.max_mdl_variants]):
                    setattr(lims_record, "drug_target_value", "Predicted low-level resistance to rifampicin. May test susceptible by phenotypic methods")
                else:
                    setattr(lims_record, "drug_target_value", "Predicted resistance to rifampicin")

                logger.debug(f"{lims_record}")
                return

            elif max_gene_code.max_mdl_interpretation in ["S"]:
                if any([self._is_synonymous_rpob_rrdr(v) for v in max_gene_code.max_mdl_variants]):
                    setattr(lims_record, "drug_target_value", "Predicted susceptibility to rifampicin. The detected synonymous mutation(s) do not confer resistance")
                else:
                    setattr(lims_record, "drug_target_value", "Predicted susceptibility to rifampicin")
                logger.debug(f"{lims_record}")
                return

        # all other drug target logic
        if max_gene_code.max_mdl_interpretation in ["R"]:
            setattr(lims_record, "drug_target_value", f"Mutation(s) associated with resistance to {lims_record.drug} detected")

        elif max_gene_code.max_mdl_interpretation in ["U"]:
            setattr(lims_record, "drug_target_value", f"The detected mutation(s) have uncertain significance. Resistance to {lims_record.drug} cannot be ruled out")

        elif max_gene_code.max_mdl_interpretation in ["S", "WT", "NA"]:
            setattr(lims_record, "drug_target_value", f"No mutations associated with resistance to {lims_record.drug} detected")

        elif max_gene_code.max_mdl_interpretation in ["Insufficient Coverage"]:
            setattr(lims_record, "drug_target_value", f"Pending Retest")

        logger.debug(f"Resolved drug target value: {lims_record}")
        return


    def set_gene_code_max_mdl(self, gene_code: LIMSGeneCode, variants: list[Variant]) -> None:
        """Set and assign the highest MDL interpretation from a list of Variants to a LIMSGeneCode object

        Args:
            variants: List of Variant objects

        Returns:
            - The highest-ranked MDL interpretation string, or "NA" if no variants have an interpretation.
            - List of Variants that have this highest MDL interpretation (could be multiple if tied)
        """
        all_mdl_interpretations = [v.mdl_interpretation for v in variants if v.mdl_interpretation]

        max_mdl_interpretation = max(
            all_mdl_interpretations,
            key=lambda x: self.RESISTANCE_RANKING.get(x, -1),
            default="NA"
        )

        max_mdl_variants = [v for v in variants if v.mdl_interpretation == max_mdl_interpretation]

        # reportable variants: R, U, and synonymous RRDR (drawn from all variants)
        max_mdl_reportable_variants = [
            v for v in variants
            if v.mdl_interpretation in ["R", "U"] or self._is_synonymous_rpob_rrdr(v)
        ]

        # assigning max_mdl_variants/interpretation to LIMSGeneCode object
        setattr(gene_code, "max_mdl_interpretation", max_mdl_interpretation)
        setattr(gene_code, "max_mdl_variants", max_mdl_variants)
        setattr(gene_code, "max_mdl_reportable_variants", max_mdl_reportable_variants)

        logger.debug(f"Setting max MDL: {gene_code}")
        return


    def _is_synonymous_rpob_rrdr(self, variant: Variant) -> bool:
        # Codon/nucleotide positions for genes requiring special LIMS handling
        SPECIAL_POSITIONS = {
            "rpoB": [426, 452],
        }

        return (
            variant.gene_name == "rpoB" and
            variant.drug == "rifampicin" and
            variant.type == "synonymous_variant" and
            Helper.is_mutation_within_range(
                Helper.get_position(variant.protein_change),
                SPECIAL_POSITIONS["rpoB"]
            )
        )

    def resolve_gene_target(
        self,
        gene_code: LIMSGeneCode,
    ) -> None:
        """Resolve gene text for S max MDL interpretation.

        For rpoB + rifampicin: checks for RRDR synonymous variants.
        For other genes: checks if gene itself has S interpretation.

        Args:
            gene_code: LIMSGeneCode object for this specific gene

        Returns:
            String of gene_target_value to be assigned for given LIMSGeneCode object
        """
        # only R, U, and synonymous RRDR variants are reportable on the gene target
        if gene_code.max_mdl_reportable_variants:
            # never report "NA" protein_changes; this is different behavior than lab report
            gene_target_value = "; ".join(
                f"{v.protein_change if v.protein_change != 'NA' else v.nucleotide_change}"
                f"{' [synonymous]' if self._is_synonymous_rpob_rrdr(v) else ''}"
                for v in gene_code.max_mdl_reportable_variants
            )
            setattr(gene_code, "gene_target_value", gene_target_value)

        elif gene_code.max_mdl_interpretation in ["S"]:
            setattr(gene_code, "gene_target_value", "No high confidence mutations detected")

        elif gene_code.max_mdl_interpretation in ["WT", "NA"]:
            setattr(gene_code, "gene_target_value", "No mutations detected")

        # overwrite mutation list if max_mdl_interpretation if Insufficient Coverage
        if gene_code.max_mdl_interpretation in ["Insufficient Coverage"]:
            setattr(gene_code, "gene_target_value", "No sequence")

        logger.debug(f"Resolved gene target value: {gene_code}")
        return

    def _passes_lims_coverage_fraction(
        self,
        lims_records: list[LIMSRecord],
        locus_coverage_map: Dict[str, LocusCoverage],
    ) -> bool:
        """Compute fraction of LIMS genes passing coverage QC.

        Args:
            lims_records: List of LIMSRecord objects to check
            locus_coverage_map: Mapping of locus_tag names to their coverage percentages

        Returns:
            bool: True if fraction of LIMS genes passing coverage QC is above threshold, False otherwise
        """
        lims_genes = set()
        # get list of all unique genes found in the LIMS records
        for rec in lims_records:
            for gene in rec.gene_codes:
                lims_genes.add(gene)

        # Only count/check coverage for the locus tag corresponding to the genes found in the LIMS format
        passing_genes = 0
        for lims_gene in lims_genes:

            # These should never happen, we should always find a locus tag; see check_bed_for_lims_genes in check_inputs.py
            lims_locus_tag = GeneDatabase.get_locus_tag(lims_gene)
            if not lims_locus_tag:
                raise ValueError(f"Could not find locus tag for LIMS gene {lims_gene} in GeneDatabase. Check that the input bed file contains the correct locus tags for LIMS genes")
            locus_coverage = locus_coverage_map.get(lims_locus_tag)
            if not locus_coverage:
                raise ValueError(f"Could not find locus tag for LIMS gene {lims_gene} in locus coverage map. Check that the input bed file contains the correct locus tags for LIMS genes")

            # if ERR coverage exists AND the `--use_err_for_qc` flag is set, use ERR coverage for determining lims QC
            if self.config.USE_ERR_FOR_QC and locus_coverage.err_coverage:
                logger.debug(f"Using ERR coverage for LIMS QC of gene {lims_gene}|{lims_locus_tag}")
                locus_coverage = locus_coverage.err_coverage

            # Determine if locus coverage is below threshold for LIMS QC fail and check if gene contains variants with valid deletions
            has_low_boc = locus_coverage.has_breadth_below(self.config.MIN_PERCENT_COVERAGE)
            boc = locus_coverage.breadth_of_coverage
            has_valid_deletion = locus_coverage.contains_loci_with_valid_deletion(lims_locus_tag)

            if has_low_boc and not has_valid_deletion:
                logger.debug(
                    f"LIMS coverage QC FAILED for gene {lims_gene}|{lims_locus_tag}: "
                    f"BOC: {(boc):.3f} BELOW threshold {self.config.MIN_PERCENT_COVERAGE:.3f}"
                )
            else:
                logger.debug(
                    f"LIMS coverage QC PASSED for gene {lims_gene}|{lims_locus_tag}: "
                    f"BOC: {(boc):.3f} ABOVE threshold {self.config.MIN_PERCENT_COVERAGE:.3f}"
                )
                passing_genes += 1

        if (passing_genes / len(lims_genes) < self.config.MIN_PERCENT_LOCI_COVERED):
            logger.debug(
                f"LIMS coverage fraction {passing_genes}/{len(lims_genes)} = {(passing_genes / len(lims_genes)):.3f} "
                f"BELOW threshold {self.config.MIN_PERCENT_LOCI_COVERED:.3f}"
            )
            return False
        else:
            logger.debug(
                f"LIMS coverage fraction {passing_genes}/{len(lims_genes)} = {(passing_genes / len(lims_genes)):.3f} "
                f"ABOVE threshold {self.config.MIN_PERCENT_LOCI_COVERED:.3f}"
            )
            return True

    def _get_pnca_his57asp_variants(self, variants: list[Variant]) -> list[Variant]:
        """Check if any pncA His57Asp variants are present in the list of variants.

        Args:
            variants: List of Variant objects
        Returns:
            list[Variant]: List of pncA His57Asp variants
        """
        pnca_variants = []
        for variant in variants:
            if variant.gene_name == "pncA" and variant.protein_change == "p.His57Asp":
                pnca_variants.append(variant)
        return pnca_variants


    # should happen last after setting all other fields
    def process_lims_mtbc_id(
        self,
        lims_records: list[LIMSRecord],
        variants: list[Variant],
        locus_coverage_map: Dict[str, LocusCoverage],
        detected_lineage: str,
        detected_sublineage: str,
    ) -> str:
        """Determine the lineage for the LIMS report.

        Args:
            lims_records: List of LIMSRecord objects
            variants: List of all Variant objects (including WT)
            locus_coverage_map: Mapping of locus_tag names to their coverage percentages
            detected_lineage: Raw lineage string from TBProfiler JSON
            detected_sublineage: Raw sublineage string from TBProfiler JSON

        Returns:
            str: detected lineage for LIMS report
        """
        lineage = set()

        pnca_his57asp_variants = self._get_pnca_his57asp_variants(variants)

        # Can't determine lineage unless sufficient coverage of LIMS genes
        if self._passes_lims_coverage_fraction(lims_records, locus_coverage_map):
            if self.config.TNGS:
                # Section 5.5: tNGS lineage based on pncA His57Asp
                if pnca_his57asp_variants:
                    if any([v.fails_positional_qc for v in pnca_his57asp_variants]):
                        lineage.add("DNA of Mycobacterium tuberculosis complex detected (M. bovis not ruled out)")
                    else:
                        lineage.add("DNA of Mycobacterium bovis detected")
                else:
                    lineage.add("DNA of Mycobacterium tuberculosis complex detected (not M. bovis)")
            else:
                # Section 5.4: WGS lineage
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

        return "; ".join(lineage)
