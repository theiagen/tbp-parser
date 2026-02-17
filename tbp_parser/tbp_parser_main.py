import logging
from arguments import parse_arguments
from utils import (
    Configuration,
    setup_logger,
    check_dependency_exists,
    check_bed_for_lims_genes,
)
from coverage import (
    CoverageCalculator,
    parse_bed_file,
)
from variant import (
    VariantProcessor,
    VariantInterpreter,
    VariantQC,
    parse_tbprofiler_json,
)
from lims import (
    LIMSProcessor,
    parse_lims_yml_file,
)
from reporters import (
    write_laboratorian_report,
    write_lims_report,
    write_looker_report,
    write_target_coverage_report,
    write_locus_coverage_report,
)

def main():
    options = parse_arguments()
    setup_logger(
        output_prefix=options.output_prefix,
        level=logging.DEBUG if options.debug else logging.INFO,
    )

    check_dependency_exists()

    config = Configuration(options)

    # Coverage calculation
    coverage_calculator = CoverageCalculator(config)
    bed_records = parse_bed_file(config.tbdb_bed)
    bed_records = coverage_calculator.populate_reads_by_position(bed_records)
    bed_records = coverage_calculator.resolve_overlapping_regions(bed_records)

    GENE_COVERAGE_MAP, LOCUS_COVERAGE_MAP = coverage_calculator.generate_coverage_maps(bed_records)
    WILDTYPE_CANDIDATES = [_ for _ in LOCUS_COVERAGE_MAP.keys()]

    # VariantRecord parsing
    variant_records, SAMPLE_ID = parse_tbprofiler_json(config.input_json)

    # Variant processing: expansion, extraction, deduplication, unreported variant generation
    variant_processor = VariantProcessor()
    all_variants = variant_processor.process_variant_records(variant_records)
    all_variants = variant_processor.deduplicate_variants(all_variants)
    unreported_variants = variant_processor.generate_unreported_variants(all_variants, SAMPLE_ID)

    # Interpretation for all_variants and unreported_variants (defining WT/NA interpretations)
    variant_interpreter = VariantInterpreter()
    all_variants = variant_interpreter.determine_interpretation(all_variants)
    #unreported_variants = variant_interpreter.determine_interpretation(unreported_variants)

    # QC for all_variants (not necessary for unreported_variants)
    variant_qc = VariantQC(config)
    all_variants = variant_qc.apply_qc(all_variants, LOCUS_COVERAGE_MAP)
    unreported_variants = variant_qc.apply_wildtype_qc(unreported_variants, LOCUS_COVERAGE_MAP)

    # merge all_varaints and unreported variants for reporting
    all_variants = all_variants + unreported_variants

    # Determine which genes have low depth of coverage
    low_depth_genes = [
        locus_tag for locus_tag, coverage in LOCUS_COVERAGE_MAP.items()
        if coverage.has_breadth_below(config.MIN_PERCENT_COVERAGE)
    ]

    GENES_WITH_VALID_DELETIONS = [str(variant.gene_name) for variant in all_variants if variant.is_valid_deletion]

    # Write reports
    write_laboratorian_report(config, all_variants)

    _, raw_lineage, lineage_english = write_lims_report(
        config, all_variants, low_depth_genes, GENES_WITH_VALID_DELETIONS
    )

    write_looker_report(
        config, all_variants, low_depth_genes, GENES_WITH_VALID_DELETIONS,
        raw_lineage, lineage_english
    )

    write_gene_coverage_report(config, SAMPLE_ID, GENE_COVERAGE_MAP)
    write_locus_coverage_report(config, SAMPLE_ID, LOCUS_COVERAGE_MAP)


if __name__ == "__main__":
    main()
