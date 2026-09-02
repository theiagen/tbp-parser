import logging

from tbp_parser.arguments import (
    BUILD_GENE_DB_COMMAND,
    BUILD_LIMS_FMT_COMMAND,
    PARSE_COMMAND,
    parse_arguments,
)
from tbp_parser.GeneDB import GeneDatabase, build_gene_db
from tbp_parser.Utilities import (
    Configuration,
    setup_logger,
    validate_inputs,
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
    build_lims_fmt,
    parse_lims_yml_file,
)
from tbp_parser.Reporters import (
    write_laboratorian_report,
    write_lims_report,
    write_looker_report,
    write_coverage_report,
)

logger = logging.getLogger(__name__)

def main():
    options = parse_arguments()

    setup_logger(
        log_file=getattr(options, "output_prefix", None),
        level=logging.DEBUG if options.debug else logging.INFO,
    )

    logger.info(f"\n\nExecuting subcommand: {options.command}\n")

    if options.command == BUILD_GENE_DB_COMMAND:
        db_bed_records = parse_bed_file(options.db_bed)
        build_gene_db(db_bed_records, options.output)
        return

    if options.command == BUILD_LIMS_FMT_COMMAND:
        build_lims_fmt(options.gene_database_yml, options.output)
        return

    if options.command == PARSE_COMMAND:
        parse(options)
        return

def parse(options):
    # Set up Configuration and GeneDatabase instances to be used throughout
    config = Configuration(options)
    gdb = GeneDatabase(config.gene_database_yml)

    # Parse input files
    lims_records = parse_lims_yml_file(config.lims_report_format_yml)
    bed_records = parse_bed_file(config.coverage_bed)
    err_records = parse_bed_file(config.err_coverage_bed)
    variant_records, SAMPLE_ID, LINEAGE_ID, SUBLINEAGE_ID = parse_tbprofiler_json(config.input_json)

    # Validate all gene/drug associations from the input files are present in the Gene Database
    if not config.SKIP_INPUT_VALIDATION:
        validate_inputs(
            bed_records=bed_records,
            lims_records=lims_records,
            variant_records=variant_records,
        )

    # Coverage calculation
    coverage_calculator = CoverageCalculator()
    LOCUS_COVERAGE_MAP, TARGET_COVERAGE_MAP = coverage_calculator.calculate(bed_records, err_records)

    # Variant processing: expansion, extraction, deduplication, unreported variant generation
    variant_processor = VariantProcessor()
    reported_variants, unreported_variants = variant_processor.process(variant_records, SAMPLE_ID)

    # Interpretation for reported_variants (not needed for unreported_variants)
    variant_interpreter = VariantInterpreter()
    reported_variants = variant_interpreter.determine_interpretation(reported_variants)

    # QC for reported_variants and unreported_variants
    variant_qc = VariantQC()
    combined_variants = variant_qc.qc(
        variants=reported_variants,
        unreported_variants=unreported_variants,
        locus_coverage_map=LOCUS_COVERAGE_MAP,
        target_coverage_map=TARGET_COVERAGE_MAP,
    )

    # Process all LIMS records and lineage information for final report
    lims_processor = LIMSProcessor()
    lims_records, lims_lineage = lims_processor.process(
        lims_records=lims_records,
        variants=combined_variants,
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
    write_laboratorian_report(combined_variants)

    # Write Looker report
    write_looker_report(
        variants=combined_variants,
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
