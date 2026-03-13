import logging
from tbp_parser.arguments import parse_arguments
from tbp_parser.Utilities import (
    Configuration,
    GeneDatabase,
    setup_logger,
    check_bed_for_lims_genes,
)
from tbp_parser.Coverage import (
    CoverageCalculator,
    parse_bed_file,
)
from tbp_parser.Variant import (
    VariantProcessor,
    VariantInterpreter,
    VariantQC,
    parse_tbprofiler_json,
)
from tbp_parser.LIMS import (
    LIMSProcessor,
    parse_lims_yml_file,
)
from tbp_parser.Reporters import (
    write_laboratorian_report,
    write_lims_report,
    write_looker_report,
    write_coverage_report,
)

def main():
    options = parse_arguments()
    setup_logger(
        output_prefix=options.output_prefix,
        level=logging.DEBUG if options.debug else logging.INFO,
    )

    # Set up Configuration and GeneDatabase instances to be used throughout
    config = Configuration(options)
    gdb = GeneDatabase(config.gene_database_yml)

    # Check entries match between LIMS and BED input files before processing
    lims_records = parse_lims_yml_file(config.lims_report_format_yml)
    bed_records = parse_bed_file(config.coverage_bed)
    err_records = parse_bed_file(config.err_coverage_bed)
    check_bed_for_lims_genes(bed_records, lims_records)

    # Coverage calculation
    coverage_calculator = CoverageCalculator()
    LOCUS_COVERAGE_MAP, TARGET_COVERAGE_MAP = coverage_calculator.calculate(bed_records, err_records)

    # VariantRecord parsing
    variant_records, SAMPLE_ID, LINEAGE_ID, SUBLINEAGE_ID = parse_tbprofiler_json(config.input_json)

    # Variant processing: expansion, extraction, deduplication, unreported variant generation
    variant_processor = VariantProcessor()
    reported_variants, unreported_variants = variant_processor.process(variant_records, SAMPLE_ID)

    # Interpretation for reported_variants (not needed for unreported_variants)
    variant_interpreter = VariantInterpreter()
    reported_variants = variant_interpreter.determine_interpretation(reported_variants)

    # QC for reported_variants and unreported_variants
    variant_qc = VariantQC()
    reported_variants = variant_qc.qc(
        variants=reported_variants,
        unreported_variants=unreported_variants,
        locus_coverage_map=LOCUS_COVERAGE_MAP,
        target_coverage_map=TARGET_COVERAGE_MAP,
    )

    # Process all LIMS records and lineage information for final report
    lims_processor = LIMSProcessor()
    lims_records, lims_lineage = lims_processor.process(
        lims_records=lims_records,
        variants=reported_variants,
        locus_coverage_map=LOCUS_COVERAGE_MAP,
        detected_lineage=LINEAGE_ID,
        detected_sublineage=SUBLINEAGE_ID,
    )

    # Write LIMS report
    write_lims_report(
        lims_records=lims_records,
        lims_lineage=lims_lineage,
        sample_id=SAMPLE_ID,
        detected_lineage=LINEAGE_ID,
    )

    # Write lab report
    write_laboratorian_report(reported_variants)

    # Write Looker report
    write_looker_report(
        variants=reported_variants,
        lims_lineage=lims_lineage,
        sample_id=SAMPLE_ID,
        detected_lineage=LINEAGE_ID,
    )

    # Write coverage reports
    write_coverage_report(
        sample_name=SAMPLE_ID,
        coverage_map=LOCUS_COVERAGE_MAP,
    )

    # Only write target coverage report if there are more targets than loci (i.e. some genes have multiple/split targets)
    # Otherwise the target coverage report would be redundant with the locus coverage report
    if len(TARGET_COVERAGE_MAP) > len(LOCUS_COVERAGE_MAP):
        write_coverage_report(
            sample_name=SAMPLE_ID,
            coverage_map=TARGET_COVERAGE_MAP,
        )

if __name__ == "__main__":
    main()
