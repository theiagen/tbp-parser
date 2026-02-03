import logging
from arguments import parse_arguments

from utils.check_inputs import check_dependency_exists
from utils.config import Configuration
from utils.logger_setup import setup_logger

from coverage import parse_bed_file, CoverageCalculator
from variant import parse_tbprofiler_json, VariantProcessor, VariantInterpreter, VariantQC

from reporters import (
    write_laboratorian_report,
    write_lims_report,
    write_looker_report,
    write_gene_coverage_report,
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

    GENE_COVERAGE_MAP = coverage_calculator.generate_gene_coverage_map(bed_records)
    LOCUS_COVERAGE_MAP = coverage_calculator.generate_locus_coverage_map(bed_records)
    WILDTYPE_CANDIDATES = [_ for _ in LOCUS_COVERAGE_MAP.keys()]

    # Variant parsing, deduplication, and unreported variant generation
    all_variants, SAMPLE_ID = parse_tbprofiler_json(config.input_json)
    all_variants = VariantProcessor.deduplicate_variants(all_variants)
    unreported_variants = VariantProcessor.generate_unreported_variants(all_variants, SAMPLE_ID, WILDTYPE_CANDIDATES)

    # Interpretation for all_variants and unreported_variants (defining WT/NA interpretations)
    variant_interpreter = VariantInterpreter()
    all_variants = variant_interpreter.determine_interpretation(all_variants)
    unreported_variants = variant_interpreter.determine_interpretation(unreported_variants)

    # QC for all_variants (not necessary for unreported_variants)
    variant_qc = VariantQC(config)
    variant_qc.add_qc_warning(all_variants, LOCUS_COVERAGE_MAP)
    genes_with_valid_deletions = VariantQC.get_genes_with_valid_deletions(all_variants)

    # Determine which genes have low depth of coverage
    low_depth_genes = [
        locus_tag for locus_tag, coverage in LOCUS_COVERAGE_MAP.items()
        if coverage.has_breadth_below(config.MIN_PERCENT_COVERAGE)
    ]

    reporter = Reporter(config)
    reporter.write_laboratorian_report(all_variants)

    # Write reports
    write_laboratorian_report(config, all_variants)

    # use err option here if applicable

    # laboratorian = Laboratorian(config, coverage.bed_records)
    #laboratorian.run()


if __name__ == "__main__":
    main()
