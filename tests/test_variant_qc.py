from variant import VariantQC, QCResult
from coverage import LocusCoverage


class TestIsDeletion:
    def test_deletion_detected(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(nucleotide_change="c.1_100del")
        assert qc._is_deletion(v) is True

    def test_non_deletion(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(nucleotide_change="c.1349C>T")
        assert qc._is_deletion(v) is False


class TestFailsTngsSpecificQc:
    def test_rrs_low_frequency_fails(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(gene_name="rrs", freq=0.05, depth=100)
        assert qc._fails_tngs_specific_qc(v) is True

    def test_rrs_low_read_support_fails(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(gene_name="rrs", freq=0.50, depth=5)
        assert qc._fails_tngs_specific_qc(v) is True

    def test_rrs_good_qc_passes(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(gene_name="rrs", freq=0.50, depth=100)
        assert qc._fails_tngs_specific_qc(v) is False

    def test_rrl_low_frequency_fails(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(gene_name="rrl", freq=0.05, depth=100)
        assert qc._fails_tngs_specific_qc(v) is True

    def test_non_special_gene_passes(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(gene_name="katG", freq=0.50, depth=100)
        assert qc._fails_tngs_specific_qc(v) is False


class TestFailsTngsBoundaryQc:
    def test_low_read_support_fails(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(freq=0.50, depth=5)  # read_support = 2.5
        assert qc._fails_tngs_boundary_qc(v) is True

    def test_mid_read_support_low_freq_fails(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(freq=0.05, depth=1000)  # RS=50 in [10,100), freq < 0.25
        assert qc._fails_tngs_boundary_qc(v) is True

    def test_high_read_support_low_freq_fails(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(freq=0.05, depth=5000)  # RS=250 >= 100, freq < 0.10
        assert qc._fails_tngs_boundary_qc(v) is True

    def test_high_read_support_good_freq_passes(self, mock_config, make_variant):
        qc = VariantQC(mock_config)
        v = make_variant(freq=0.50, depth=1000)  # RS=500 >= 100, freq >= 0.10
        assert qc._fails_tngs_boundary_qc(v) is False


class TestApplyQc:
    def test_low_depth_fails_positional_qc(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(depth=5, freq=0.50)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is True
        assert "Failed quality in the mutation position" in v.warning

    def test_low_frequency_fails_positional_qc(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(depth=100, freq=0.01)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is True
        assert "Failed quality in the mutation position" in v.warning

    def test_good_variant_passes_qc(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(depth=100, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is not True
        assert len(v.warning) == 0

    def test_low_breadth_adds_locus_warning_for_non_r(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(depth=100, freq=0.95, confidence="Uncertain significance")
        # Pre-set interpretation as determine_interpretation would
        v.mdl_interpretation = "U"
        v.looker_interpretation = "U"
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=50.0)
        qc.apply_qc([v], {"Rv0667": locus})

        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"

    def test_low_breadth_r_mutation_no_positional_fail_keeps_interpretation(
        self, mock_config, make_variant, make_locus_coverage
    ):
        qc = VariantQC(mock_config)
        v = make_variant(depth=100, freq=0.95, confidence="Assoc w R")
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=50.0)
        qc.apply_qc([v], {"Rv0667": locus})

        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "R"  # preserved for R without positional fail

    def test_tngs_outside_region(self, mock_config, make_variant, make_locus_coverage):
        mock_config.TNGS = True
        qc = VariantQC(mock_config)
        v = make_variant(depth=100, freq=0.95, pos=1)  # pos outside locus coords
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(759807, 763325)])
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is True
        assert "This mutation is outside the expected region" in v.warning
        assert v.mdl_interpretation == "NA"

    def test_deletion_with_good_qc_is_valid(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.is_valid_deletion is True
        assert v.fails_qc is not True

    def test_deletion_with_zero_depth_good_freq_passes(self, mock_config, make_variant, make_locus_coverage):
        """Rule 4.2.1.3: depth=0 with good frequency should pass (TB Profiler quirk)."""
        qc = VariantQC(mock_config)
        v = make_variant(nucleotide_change="c.1_100del", depth=0, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is not True
        assert v.is_valid_deletion is True

    def test_deletion_with_low_depth_fails(self, mock_config, make_variant, make_locus_coverage):
        """Rule 4.2.1.2: deletion with some depth but below threshold fails."""
        qc = VariantQC(mock_config)
        v = make_variant(nucleotide_change="c.1_100del", depth=5, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is True
        assert "Failed quality in the mutation position" in v.warning
        assert v.is_valid_deletion is False

    def test_r_mutation_both_positional_and_locus_fail_overwrite(
        self, mock_config, make_variant, make_locus_coverage
    ):
        """Rule 4.2.2.3.4: R mutation with both positional and locus fail overwrites."""
        qc = VariantQC(mock_config)
        v = make_variant(depth=5, freq=0.50, confidence="Assoc w R")
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=50.0)
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.fails_qc is True
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"
        assert v.looker_interpretation == "Insufficient Coverage"


class TestGetGenesWithValidDeletions:
    def test_valid_deletion_tracked(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(
            nucleotide_change="c.1_100del",
            depth=100, freq=0.95,
        )
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.is_valid_deletion is True

    def test_failed_deletion_not_tracked(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(
            nucleotide_change="c.1_100del",
            depth=2, freq=0.01,  # fails QC
        )
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.is_valid_deletion is False

    def test_non_deletion_not_tracked(self, mock_config, make_variant, make_locus_coverage):
        qc = VariantQC(mock_config)
        v = make_variant(depth=100, freq=0.95)  # SNP, not deletion
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus})

        assert v.is_valid_deletion is False
