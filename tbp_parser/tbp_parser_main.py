import logging
from arguments import parse_arguments

from utils.check_inputs import check_dependency_exists
from utils.config import Configuration
from utils.logger_setup import setup_logger


from coverage import parse_bed_file
from variant import parse_tbprofiler_json

from coverage import CoverageCalculator
from variant import VariantProcessor
from variant import VariantInterpreter
from reporters.reporter import Reporter



def main():
    options = parse_arguments()
    setup_logger(logging.DEBUG if options.debug else logging.INFO)

    check_dependency_exists()

    config = Configuration(options)

    coverage_calculator = CoverageCalculator(config)
    bed_records = parse_bed_file(config.tbdb_bed)
    bed_records = coverage_calculator.populate_reads_by_position(bed_records)
    bed_records = coverage_calculator.resolve_overlapping_regions(bed_records)

    GENE_COVERAGE_MAP = coverage_calculator.generate_gene_coverage_map(bed_records)
    LOCUS_COVERAGE_MAP = coverage_calculator.generate_locus_coverage_map(bed_records)

    all_variants = parse_tbprofiler_json(config.input_json)
    all_variants = VariantProcessor.deduplicate_variants(all_variants)
    all_variants += VariantProcessor.generate_unreported_variants(all_variants)

    VariantInterpreter.determine_interpretation(all_variants)
    # VariantQC.perform_variant_qc(all_variants, GENE_COVERAGE_LIST, LOCUS_COVERAGE_LIST)

    reporter = Reporter(config)
    reporter.write_laboratorian_report(all_variants)

    breakpoint()

    # use err option here if applicable

    # laboratorian = Laboratorian(config, coverage.bed_records)
    #laboratorian.run()


if __name__ == "__main__":
    main()