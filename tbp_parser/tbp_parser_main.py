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

    config = Configuration(options)

    # Perform necessary input checks before processing
    check_dependency_exists()

    # Check entries match between LIMS and BED input files before processing
    lims_records = parse_lims_yml_file(config.lims_report_format_yml)
    bed_records = parse_bed_file(config.tbdb_bed)
    check_bed_for_lims_genes(bed_records, lims_records)

    # Coverage calculation
    coverage_calculator = CoverageCalculator(config)
    TARGET_COVERAGE_MAP, LOCUS_COVERAGE_MAP = coverage_calculator.calculate(bed_records)

    # VariantRecord parsing
    variant_records, SAMPLE_ID, LINEAGE_ID, SUBLINEAGE_ID = parse_tbprofiler_json(config.input_json)

    # Variant processing: expansion, extraction, deduplication, unreported variant generation
    variant_processor = VariantProcessor()
    all_variants, unreported_variants = variant_processor.process(variant_records, SAMPLE_ID)

    # Interpretation for all_variants (not needed for unreported variants)
    variant_interpreter = VariantInterpreter()
    all_variants = variant_interpreter.determine_interpretation(all_variants)

    # QC for all_variants and unreported variants
    variant_qc = VariantQC(config)
    all_variants, GENES_WITH_VALID_DELETIONS = variant_qc.qc(
        variants=all_variants,
        unreported_variants=unreported_variants,
        locus_coverage_map=LOCUS_COVERAGE_MAP,
    )

    # Write lab report
    write_laboratorian_report(config, all_variants)

    # Process all LIMS records and lineage information for final report
    lims_processor = LIMSProcessor(config)
    lims_records = lims_processor.process_lims_records(lims_records, all_variants)
    lims_lineage = lims_processor.process_lims_mtbc_id(
        lims_records=lims_records,
        variants=all_variants,
        locus_coverage_map=LOCUS_COVERAGE_MAP,
        genes_with_valid_deletions=GENES_WITH_VALID_DELETIONS,
        detected_lineage=LINEAGE_ID,
        detected_sublineage=SUBLINEAGE_ID,
    )

    # Write LIMS report
    write_lims_report(
        config=config,
        lims_records=lims_records,
        lims_lineage=lims_lineage,
        sample_id=SAMPLE_ID,
        detected_lineage=LINEAGE_ID,
    )

    # Write Looker report
    write_looker_report(
        config=config,
        variants=all_variants,
        lims_lineage=lims_lineage,
        sample_id=SAMPLE_ID,
        detected_lineage=LINEAGE_ID,
    )

    # Write coverage reports
    write_target_coverage_report(
        config,
        sample_name=SAMPLE_ID,
        target_coverage_map=TARGET_COVERAGE_MAP,
        genes_with_valid_deletions=GENES_WITH_VALID_DELETIONS,
    )
    write_locus_coverage_report(
        config,
        sample_name=SAMPLE_ID,
        locus_coverage_map=LOCUS_COVERAGE_MAP,
        genes_with_valid_deletions=GENES_WITH_VALID_DELETIONS,
    )


if __name__ == "__main__":
    main()
