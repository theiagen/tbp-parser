import pytest
from unittest.mock import MagicMock
from variant import Variant
from coverage import LocusCoverage, GeneCoverage


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
        )
        return v
    return _make


@pytest.fixture
def make_wt_variant():
    """Factory fixture to create WT (wild-type) Variant objects."""
    def _make(gene_name="rpoB", gene_id="Rv0667", drug="rifampicin", sample_id="test_sample"):
        v = Variant.from_thin_air(sample_id=sample_id, gene_id=gene_id, gene_name=gene_name, drug=drug)
        v.mdl_interpretation = "WT"
        v.looker_interpretation = "WT"
        v.rationale = "NA"
        return v
    return _make


@pytest.fixture
def make_locus_coverage():
    """Factory fixture to create LocusCoverage objects."""
    def _make(
        locus_tag="Rv0667",
        gene_names=None,
        coords=None,
        breadth_of_coverage=100.0,
        average_depth=200.0,
    ):
        return LocusCoverage(
            locus_tag=locus_tag,
            gene_names=gene_names or ["rpoB"],
            coords=coords or [(759807, 763325)],
            breadth_of_coverage=breadth_of_coverage,
            average_depth=average_depth,
        )
    return _make


@pytest.fixture
def make_gene_coverage():
    """Factory fixture to create GeneCoverage objects."""
    def _make(
        locus_tag="Rv0667",
        gene_name="rpoB",
        coords=(759807, 763325),
        breadth_of_coverage=100.0,
        average_depth=200.0,
    ):
        return GeneCoverage(
            locus_tag=locus_tag,
            gene_name=gene_name,
            coords=coords,
            breadth_of_coverage=breadth_of_coverage,
            average_depth=average_depth,
        )
    return _make


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
    config.TNGS_SPECIFIC_QC_OPTIONS = {
        "RRS_FREQUENCY": 0.25,
        "RRS_READ_SUPPORT": 10,
        "RRL_FREQUENCY": 0.25,
        "RRL_READ_SUPPORT": 10,
        "ETHA237_FREQUENCY": 0.10,
        "RPOB449_FREQUENCY": 0.10,
    }
    return config
