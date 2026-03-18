import pytest
import pysam
from pathlib import Path
from copy import deepcopy
from unittest.mock import MagicMock
from tbp_parser.Utilities import Configuration, GeneDatabase
from tbp_parser.Variant import Variant, VariantRecord, Annotation, Consequences
from tbp_parser.Coverage import LocusCoverage, TargetCoverage, BedRecord, CoverageCalculator

@pytest.fixture
def mock_config():
    """Create a mock Configuration object with sensible defaults."""
    config = MagicMock()
    config.MIN_DEPTH = 10
    config.MIN_FREQUENCY = 0.10
    config.MIN_READ_SUPPORT = 10
    config.MIN_PERCENT_COVERAGE = 0.90
    config.MIN_PERCENT_LOCI_COVERED = 0.50
    config.TNGS = False
    config.DO_NOT_TREAT_R_MUTATIONS_DIFFERENTLY = False
    config.SEQUENCING_METHOD = "WGS"
    config.OUTPUT_PREFIX = "/tmp/test_output"
    config.OPERATOR = "test_operator"
    config.TNGS_READ_SUPPORT_BOUNDARIES = [10, 100]
    config.TNGS_FREQUENCY_BOUNDARIES = [0.10, 0.25]
    config.USE_ERR_AS_BRR = False
    config.lims_report_format_yml = str(Path(__file__).parent.parent / "tbp_parser" / "data" / "default-lims-report-format.yml")
    return config

@pytest.fixture(autouse=True)
def setup_config(mock_config):
    Configuration._instance = mock_config

@pytest.fixture(autouse=True, scope="session")
def setup_gene_database():
    GeneDatabase(db_path=str(Path(__file__).parent.parent / "tbp_parser" / "data" / "default-gene-database_2026-03-03.yml"))

@pytest.fixture
def make_bed_record():
    """Module-level helper to construct BedRecord objects with sensible defaults."""
    def _make(**kwargs):
        params = {"chrom": "Chromosome", "start": 100, "end": 200, "locus_tag": "Rv0000", "gene_name": "geneA"}
        params.update(kwargs)
        bed_record = BedRecord(**params)  # type: ignore
        return bed_record
    return _make


@pytest.fixture
def make_variant():
    """Factory fixture to create Variant objects with sensible defaults."""
    def _make(
        gene_name="rpoB",
        gene_id="Rv0667",
        nucleotide_change="c.1349C>T",
        protein_change="p.Ser450Leu",
        drug="rifampicin",
        confidence="Assoc w R",
        type="missense_variant",
        pos=761155,
        depth=100,
        freq=0.95,
        source="WHO",
        comment="",
        sample_id="test_sample",
        **kwargs,
    ):
        v = Variant(
            sample_id=sample_id,
            pos=pos,
            depth=depth,
            freq=freq,
            gene_id=gene_id,
            gene_name=gene_name,
            type=type,
            nucleotide_change=nucleotide_change,
            protein_change=protein_change,
            confidence=confidence,
            drug=drug,
            source=source,
            comment=comment,
            **kwargs,
        )
        return v
    return _make


@pytest.fixture
def make_locus_coverage():
    """Factory fixture to create LocusCoverage objects."""
    def _make(
        locus_tag="Rv0000",
        gene_names=["geneA", "geneB"],
        coords=[(100, 150), (250, 350)],
        breadth_of_coverage=1.0,
        average_depth=100.0,
        err_coverage=None,
        **kwargs,
    ):
        return LocusCoverage(
            locus_tag=locus_tag,
            gene_names=gene_names,
            coords=coords,
            breadth_of_coverage=breadth_of_coverage,
            average_depth=average_depth,
            err_coverage=err_coverage,
            **kwargs,
        )
    return _make


@pytest.fixture
def make_target_coverage():
    """Factory fixture to create TargetCoverage objects."""
    def _make(
        locus_tag="Rv0000",
        gene_name="geneA",
        coords=[(100, 150)],
        breadth_of_coverage=1.0,
        average_depth=100.0,
        err_coverage=None,
        **kwargs,
    ):
        return TargetCoverage(
            locus_tag=locus_tag,
            gene_name=gene_name,
            coords=coords,
            breadth_of_coverage=breadth_of_coverage,
            average_depth=average_depth,
            err_coverage=err_coverage,
            **kwargs,
        )
    return _make




@pytest.fixture
def make_annotation():
    """Factory fixture to create Annotation objects with sensible defaults."""
    def _make(
        drug="rifampin",
        confidence="No WHO annotation",
        source="",
        comment="",
    ):
        return Annotation(drug=drug, confidence=confidence, source=source, comment=comment)
    return _make


@pytest.fixture
def make_consequences(make_annotation):
    """Factory fixture to create Consequences objects with sensible defaults."""
    def _make(
        gene_id="Rv0678",
        gene_name="mmpR5",
        type="missense_variant",
        nucleotide_change="c.200C>T",
        protein_change="p.Val67Ile",
        annotation=[],
    ):
        return Consequences(
            gene_id=gene_id,
            gene_name=gene_name,
            type=type,
            nucleotide_change=nucleotide_change,
            protein_change=protein_change,
            annotation=annotation,
        )
    return _make


@pytest.fixture
def make_variant_record(make_annotation):
    """Factory fixture to create VariantRecord objects with sensible defaults."""
    def _make(
        sample_id="test_sample",
        pos=200,
        depth=100,
        freq=0.95,
        gene_id="Rv0005",
        gene_name="gyrB",
        type="missense_variant",
        nucleotide_change="c.100C>T",
        protein_change="p.Ser34Leu",
        annotation=None,
        consequences=None,
        gene_associated_drugs=None,
    ):
        if annotation is None:
            annotation = [make_annotation()]
        if consequences is None:
            consequences = []
        if gene_associated_drugs is None:
            gene_associated_drugs = [annotation.drug for annotation in annotation]
        return VariantRecord(
            sample_id=sample_id,
            pos=pos,
            depth=depth,
            freq=freq,
            gene_id=gene_id,
            gene_name=gene_name,
            type=type,
            nucleotide_change=nucleotide_change,
            protein_change=protein_change,
            annotation=annotation,
            consequences=consequences,
            gene_associated_drugs=gene_associated_drugs,
        )
    return _make


@pytest.fixture
def make_bam_file():
    """Factory fixture to create a BAM file for testing."""
    def _make(
        cov_start=0,
        cov_end=100,
        read_length=10,
    ):
        bam_file = str(Path(__file__).parent / "test.bam")
        # random 100 bases that I will repeat in this order to use as a reference.
        # reference base [0:10] 'GACAAGGACA' will be the same as [100:110] 'GACAAGGACA', etc
        ref_seq = "GACAAGGACATGACGTACGCGGCCCCGCTGTTCGTCACGGCCGAGTTCATCAACAACAACACCGGTGAGATCAAGAGCCAGACGGTGTTCATGGGATCGG"
        header = pysam.AlignmentHeader.from_dict({
            'HD': {'VN': '1.5', 'SO': 'coordinate'},
            'SQ': [{'SN': 'Chromosome', 'LN': int(cov_end - cov_start)}]
        })
        reads = []
        read_counter = 1

        for pos in range(cov_start, cov_end):
            read_name = f"read{read_counter}"
            read_counter += 1

            # R1
            r1 = pysam.AlignedSegment()
            r1.query_name = read_name
            r1.query_sequence = ''.join(ref_seq[(pos + i) % len(ref_seq)] for i in range(read_length))
            r1.flag = 99   # paired, proper pair, mate reverse strand, first in pair
            r1.reference_id = 0
            r1.reference_start = pos
            r1.mapping_quality = 60
            r1.cigartuples = [(0, read_length)]
            r1.next_reference_id = 0
            r1.next_reference_start = pos
            r1.template_length = read_length
            r1.query_qualities = pysam.qualitystring_to_array("?" * read_length)
            reads.append(r1)

            # R2 (mate)
            r2 = deepcopy(r1)
            r2.flag = 147  # paired, proper pair, self reverse strand, second in pair
            r2.template_length = -read_length
            reads.append(r2)

        reads.sort(key=lambda r: r.reference_start)
        with pysam.AlignmentFile(bam_file, "wb", header=header) as bam:
            for r in reads:
                bam.write(r)

        pysam.sort("-o", bam_file, bam_file)
        pysam.index(bam_file)
        return bam_file
    return _make


@pytest.fixture
def make_cov_calc(mock_config, make_bam_file):
    """Factory fixture to create a CoverageCalculator instance with a mock config."""
    def _make(cov_start=0, cov_end=100, read_length=10):
        mock_config.input_bam = make_bam_file(cov_start=cov_start, cov_end=cov_end, read_length=read_length) # bams are 0-based
        return CoverageCalculator()
    return _make
