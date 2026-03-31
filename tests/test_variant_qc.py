from tbp_parser.Variant import VariantQC
from tbp_parser.Coverage.coverage_data import ERRCoverage


class TestIsDeletion:
    def test_deletion_detected(self, make_variant):
        v = make_variant(nucleotide_change="c.1_100del")
        assert v._is_deletion_in_orf() is True

    def test_non_deletion(self, make_variant):
        v = make_variant(nucleotide_change="c.1349C>T")
        assert v._is_deletion_in_orf() is False


class TestFailsTngsBoundaryQc:
    def test_low_read_support_fails(self, make_variant):
        qc = VariantQC()
        v = make_variant(freq=0.50, depth=5)  # read_support = 2.5
        assert qc._fails_tngs_boundary_qc(v) is True

    def test_mid_read_support_low_freq_fails(self, make_variant):
        qc = VariantQC()
        v = make_variant(freq=0.05, depth=1000)  # RS=50 in [10,100), freq < 0.25
        assert qc._fails_tngs_boundary_qc(v) is True

    def test_high_read_support_low_freq_fails(self, make_variant):
        qc = VariantQC()
        v = make_variant(freq=0.05, depth=5000)  # RS=250 >= 100, freq < 0.10
        assert qc._fails_tngs_boundary_qc(v) is True

    def test_high_read_support_good_freq_passes(self, make_variant):
        qc = VariantQC()
        v = make_variant(freq=0.50, depth=1000)  # RS=500 >= 100, freq >= 0.10
        assert qc._fails_tngs_boundary_qc(v) is False


class TestApplyQc:
    def test_low_depth_fails_positional_qc(self, make_variant, make_locus_coverage):
        qc = VariantQC()
        v = make_variant(depth=5, freq=0.50)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is False
        assert "Failed quality in the mutation position" in v.warning

    def test_low_frequency_fails_positional_qc(self, make_variant, make_locus_coverage):
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.01)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is False
        assert "Failed quality in the mutation position" in v.warning

    def test_low_read_support_fails_positional_qc(self, make_variant, make_locus_coverage):
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, read_support=5)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is False
        assert "Failed quality in the mutation position" in v.warning

    def test_good_variant_passes_qc(self, make_variant, make_locus_coverage):
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert len(v.warning) == 0

    def test_low_breadth_u_mutation_overwrites_to_insufficient_coverage(self, make_variant, make_locus_coverage):
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Uncertain significance")
        # Pre-set interpretation as determine_interpretation would
        v.mdl_interpretation = "U"
        v.looker_interpretation = "U"
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.50)
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"
        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is True

    def test_low_breadth_r_mutation_no_positional_fail_keeps_interpretation(
        self, make_variant, make_locus_coverage
    ):
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Assoc w R")
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.50)
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "R"  # preserved for R without positional fail
        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False

    def test_tngs_outside_region(self, mock_config, make_variant, make_locus_coverage):
        mock_config.TNGS = True
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, pos=1)  # pos outside locus coords
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(759807, 763325)])
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is True
        assert "This mutation is outside the expected region" in v.warning
        assert v.mdl_interpretation == "NA"

    def test_deletion_with_good_depth_freq_passes_positional_qc(self, make_variant, make_locus_coverage):
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v._is_deletion_in_orf() is True
        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False

    def test_deletion_with_zero_depth_good_freq_passes(self, make_variant, make_locus_coverage):
        """Rule 4.2.1.3: depth=0 with good frequency should pass (TB Profiler quirk)."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=0, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert v._is_deletion_in_orf() is True

    def test_deletion_with_zero_depth_low_freq_fails(self, make_variant, make_locus_coverage):
        """Rule 4.2.1.4: depth=0 with low frequency should fail."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=0, freq=0.01)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is False
        assert "Failed quality in the mutation position" in v.warning

    def test_deletion_with_low_depth_fails(self, make_variant, make_locus_coverage):
        """Rule 4.2.1.2: deletion with some depth but below threshold fails."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=5, freq=0.95)
        locus = make_locus_coverage(locus_tag="Rv0667")
        qc.apply_qc([v], {"Rv0667": locus}, {})
        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is False
        assert "Failed quality in the mutation position" in v.warning

    def test_r_mutation_both_positional_and_locus_fail_overwrite(
        self, make_variant, make_locus_coverage
    ):
        """Rule 4.2.2.3.4: R mutation with both positional and locus fail overwrites."""
        qc = VariantQC()
        v = make_variant(depth=5, freq=0.50, confidence="Assoc w R")
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.50)
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_locus_qc is True
        assert v.fails_positional_qc is True
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"
        assert v.looker_interpretation == "Insufficient Coverage"

    def test_deletion_failing_positional_qc_not_in_valid_deletions(self, make_variant, make_locus_coverage):
        """Deletion that fails positional QC should not be assigned to valid_deletions."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=5, freq=0.95, pos=120)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)])
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is False
        assert locus.contains_variant_with_valid_deletion(v) is False

    def test_deletion_failing_qc_low_breadth_r_does_not_pass_via_valid_deletion(
        self, make_variant, make_locus_coverage
    ):
        """Deletion failing positional QC + low breadth + R -> Rule 4.2.2.3.4, not 4.2.2.2."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=5, freq=0.95, pos=120, confidence="Assoc w R")
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], breadth_of_coverage=0.50)
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is True
        assert locus.contains_variant_with_valid_deletion(v) is False
        assert v.mdl_interpretation == "Insufficient Coverage"
        assert v.looker_interpretation == "Insufficient Coverage"

    def test_deletion_failing_qc_low_breadth_s_does_not_pass_via_valid_deletion(
        self, make_variant, make_locus_coverage
    ):
        """Deletion failing positional QC + low breadth + S -> Rule 4.2.2.3.2."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=5, freq=0.95, pos=120, confidence="Uncertain significance")
        v.mdl_interpretation = "S"
        v.looker_interpretation = "S"
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], breadth_of_coverage=0.50)
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is True
        assert locus.contains_variant_with_valid_deletion(v) is False
        assert v.mdl_interpretation == "Insufficient Coverage"
        assert v.looker_interpretation == "Insufficient Coverage"

    def test_multiple_deletions_only_passing_qc_assigned(self, make_variant, make_locus_coverage):
        """Only the deletion passing positional QC should be in valid_deletions."""
        qc = VariantQC()
        good = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120, drug="rifampicin")
        bad = make_variant(nucleotide_change="c.1_100del", depth=5, freq=0.95, pos=130, drug="isoniazid")
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)])
        qc.apply_qc([good, bad], {"Rv0667": locus}, {})

        assert locus.contains_variant_with_valid_deletion(good) is True
        assert locus.contains_variant_with_valid_deletion(bad) is False
        assert good.fails_positional_qc is False
        assert good.fails_locus_qc is False
        assert bad.fails_positional_qc is True
        assert bad.fails_locus_qc is False

    def test_valid_deletion_low_breadth_passes_locus_qc(self, make_variant, make_locus_coverage):
        """Rule 4.2.2.2: Good deletion + low breadth -> passes locus QC."""
        qc = VariantQC()
        v = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], breadth_of_coverage=0.50)
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert locus.contains_variant_with_valid_deletion(v) is True

    def test_mutation_passes_locus_qc_when_locus_has_valid_deletion(self, make_variant, make_locus_coverage):
        """Rule 4.2.2.2: SNP at a locus with a sibling valid deletion + low breadth -> passes locus QC."""
        qc = VariantQC()
        # Deletion processed first so it gets assigned to valid_deletions before the SNP's locus QC
        deletion = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)
        snp = make_variant(nucleotide_change="c.1349C>T", depth=100, freq=0.95, pos=130)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], breadth_of_coverage=0.50)
        qc.apply_qc(variants=[deletion, snp], locus_coverage_map={"Rv0667": locus}, target_coverage_map={})

        assert snp.fails_positional_qc is False
        assert snp.fails_locus_qc is False
        assert locus.contains_loci_with_valid_deletion(snp.gene_id) is True
        assert locus.contains_variant_with_valid_deletion(snp) is False  # SNP itself is NOT a valid deletion

    def test_mutation_passes_locus_qc_when_locus_has_valid_deletion_but_fails_positional_qc(
        self, make_variant, make_locus_coverage
    ):
        """Rule 4.2.2.2: SNP failing positional QC at a locus with valid deletion + low breadth -> fails locus QC."""
        qc = VariantQC()
        deletion = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)
        snp = make_variant(nucleotide_change="c.1349C>T", depth=5, freq=0.50, pos=130)  # low depth -> fails positional
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], breadth_of_coverage=0.50)
        qc.apply_qc(variants=[deletion, snp], locus_coverage_map={"Rv0667": locus}, target_coverage_map={})

        assert snp.fails_positional_qc is True
        assert snp.fails_locus_qc is True
        assert "Insufficient coverage in locus" in snp.warning
        assert locus.contains_loci_with_valid_deletion(snp.gene_id) is True

    def test_wildtype_fails_locus_qc(self, make_variant, make_locus_coverage):
        """Rule 4.2.2.1: WT/NA variant should fail locus QC if breadth is low."""
        qc = VariantQC()
        v = make_variant(depth=10, freq=1, confidence="NA", nucleotide_change="NA", protein_change="NA", type="NA", pos=1)
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.50)
        qc.apply_wildtype_qc([v], {"Rv0667": locus})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is True
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"

    def test_wildtype_passes_locus_qc(self, make_variant, make_locus_coverage):
        """Rule 4.2.2.1: WT/NA variant should pass locus QC if breadth is sufficient."""
        qc = VariantQC()
        v = make_variant(depth=10, freq=1, confidence="NA", nucleotide_change="NA", protein_change="NA", type="NA", pos=1)
        locus = make_locus_coverage(locus_tag="Rv0667", breadth_of_coverage=0.95)
        qc.apply_wildtype_qc([v], {"Rv0667": locus})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert "Insufficient coverage in locus" not in v.warning
        assert v.mdl_interpretation == "WT"

class TestLocusQcWithErrCoverage:
    """Tests for locus QC when USE_ERR_FOR_QC flag swaps coverage to ERR."""

    def _make_err(self, breadth=0.95, coords=None):
        if coords is None:
            coords = [(100, 150), (250, 350)]
        return ERRCoverage(coords=coords, breadth_of_coverage=breadth, average_depth=100.0)

    def test_err_high_boc_passes_when_locus_low(self, mock_config, make_variant, make_locus_coverage):
        mock_config.USE_ERR_FOR_QC = True
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Uncertain significance")
        v.mdl_interpretation = "S"
        v.looker_interpretation = "S"
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.50,
            err_coverage=self._make_err(breadth=0.95),
        )
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert "Insufficient coverage in locus" not in v.warning
        assert v.mdl_interpretation == "S"  # preserved — ERR breadth is good

    def test_err_also_low_fails_insufficient(self, mock_config, make_variant, make_locus_coverage):
        mock_config.USE_ERR_FOR_QC = True
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Uncertain significance")
        v.mdl_interpretation = "S"
        v.looker_interpretation = "S"
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.50,
            err_coverage=self._make_err(breadth=0.50),
        )
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is True
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"

    def test_err_flag_off_uses_locus(self, mock_config, make_variant, make_locus_coverage):
        mock_config.USE_ERR_FOR_QC = False
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Uncertain significance")
        v.mdl_interpretation = "S"
        v.looker_interpretation = "S"
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.50,
            err_coverage=self._make_err(breadth=0.95),
        )
        qc.apply_qc([v], {"Rv0667": locus}, {})

        # Flag off -> uses locus breadth (0.50) -> fails
        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is True
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"

    def test_err_no_err_coverage_uses_locus(self, mock_config, make_variant, make_locus_coverage):
        mock_config.USE_ERR_FOR_QC = True
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Uncertain significance")
        v.mdl_interpretation = "S"
        v.looker_interpretation = "S"
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.95,
            err_coverage=None,  # no ERR coverage
        )
        qc.apply_qc([v], {"Rv0667": locus}, {})

        # Falls through to locus breadth (0.95) -> passes
        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert "Insufficient coverage in locus" not in v.warning
        assert v.mdl_interpretation == "S"

    def test_r_locus_fail_only_warning_preserved(self, mock_config, make_variant, make_locus_coverage):
        mock_config.USE_ERR_FOR_QC = True
        qc = VariantQC()
        v = make_variant(depth=100, freq=0.95, confidence="Assoc w R")
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.50,
            err_coverage=self._make_err(breadth=0.50),
        )
        qc.apply_qc([v], {"Rv0667": locus}, {})

        # R + locus fail only (no positional fail) -> warning added but interpretation preserved
        assert v.fails_positional_qc is False
        assert v.fails_locus_qc is False
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "R"

    def test_r_both_fail_overwrites(self, mock_config, make_variant, make_locus_coverage):
        mock_config.USE_ERR_FOR_QC = True
        qc = VariantQC()
        v = make_variant(depth=5, freq=0.50, confidence="Assoc w R")  # low depth -> positional fail
        v.mdl_interpretation = "R"
        v.looker_interpretation = "R"
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.50,
            err_coverage=self._make_err(breadth=0.50),
        )
        qc.apply_qc([v], {"Rv0667": locus}, {})

        assert v.fails_positional_qc is True
        assert v.fails_locus_qc is True
        assert "Insufficient coverage in locus" in v.warning
        assert v.mdl_interpretation == "Insufficient Coverage"
        assert v.looker_interpretation == "Insufficient Coverage"

    def test_err_valid_deletion_passes_locus_qc(self, mock_config, make_variant, make_locus_coverage, make_target_coverage):
        """Rule 4.2.2.2: Deletion in ERR valid_deletions + low ERR breadth -> passes locus QC."""
        mock_config.USE_ERR_FOR_QC = True
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)

        err = self._make_err(breadth=0.50, coords=[(100, 200)])
        locus = make_locus_coverage(
            locus_tag="Rv0667", breadth_of_coverage=0.50,
            coords=[(100, 200)],
            err_coverage=err,
        )
        target = make_target_coverage(
            locus_tag="Rv0667", gene_name="rpoB",
            coords=[(100, 200)],
            err_coverage=ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.50, average_depth=100.0),
        )

        # First assign valid deletions, then apply QC
        qc.apply_qc([del_variant], {"Rv0667": locus}, {"rpoB": target})

        assert del_variant.fails_positional_qc is False
        assert del_variant.fails_locus_qc is False


class TestAssignValidDeletions:
    """Tests for assign_variants_with_valid_deletions with and without ERR coverage."""

    def test_deletion_assigned(self, make_variant, make_locus_coverage, make_target_coverage):
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)

        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], err_coverage=None)
        target = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB", coords=[(100, 200)], err_coverage=None)

        qc.assign_variants_with_valid_deletions(del_variant, {"Rv0667": locus}, {"rpoB": target})

        assert locus.contains_variant_with_valid_deletion(del_variant) is True
        assert target.contains_variant_with_valid_deletion(del_variant) is True

    def test_deletion_outside_coverage_coords_not_assigned(self, make_variant, make_locus_coverage, make_target_coverage):
        """Good deletion at a position outside coverage coords should not be assigned."""
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=500)

        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)])
        target = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB", coords=[(100, 200)])

        qc.assign_variants_with_valid_deletions(del_variant, {"Rv0667": locus}, {"rpoB": target})

        assert locus.contains_variant_with_valid_deletion(del_variant) is False
        assert target.contains_variant_with_valid_deletion(del_variant) is False

    def test_deletion_overlapping_coverage_coords_assigned(self, make_variant, make_locus_coverage, make_target_coverage):
        """Good deletion overlapping multiple target coverage objects should be assigned to both."""
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=150)

        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)])
        target1 = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB1", coords=[(100, 200)])
        target2 = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB2", coords=[(110, 210)])

        qc.assign_variants_with_valid_deletions(del_variant, {"Rv0667": locus}, {"rpoB1": target1, "rpoB2": target2})

        assert locus.contains_variant_with_valid_deletion(del_variant) is True
        assert target1.contains_variant_with_valid_deletion(del_variant) is True
        assert target2.contains_variant_with_valid_deletion(del_variant) is True
        assert target1.valid_deletions == target2.valid_deletions == [del_variant]

    def test_non_deletion_not_assigned(self, make_variant, make_locus_coverage, make_target_coverage):
        qc = VariantQC()
        snp_variant = make_variant(depth=100, freq=0.95, pos=120)  # SNP, not deletion

        err = ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.95, average_depth=100.0)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], err_coverage=err)
        target = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB", coords=[(100, 200)])

        qc.assign_variants_with_valid_deletions(snp_variant, {"Rv0667": locus}, {"rpoB": target})

        assert locus.contains_variant_with_valid_deletion(snp_variant) is False
        assert target.contains_variant_with_valid_deletion(snp_variant) is False
        assert err.contains_variant_with_valid_deletion(snp_variant) is False

    def test_deletion_failing_positional_qc_not_assigned(self, make_variant, make_locus_coverage, make_target_coverage):
        """Deletion with fails_positional_qc=True should be skipped by assign_variants_with_valid_deletions."""
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)
        del_variant.fails_positional_qc = True

        err = ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.95, average_depth=100.0)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], err_coverage=err)
        target = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB", coords=[(100, 200)])

        qc.assign_variants_with_valid_deletions(del_variant, {"Rv0667": locus}, {"rpoB": target})

        assert locus.contains_variant_with_valid_deletion(del_variant) is False
        assert target.contains_variant_with_valid_deletion(del_variant) is False
        assert err.contains_variant_with_valid_deletion(del_variant) is False

    def test_deletion_in_err_range_assigned(self, make_variant, make_locus_coverage, make_target_coverage):
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=120)

        err_locus = ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.95, average_depth=100.0)
        err_target = ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.95, average_depth=100.0)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200)], err_coverage=err_locus)
        target = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB", coords=[(100, 200)], err_coverage=err_target)

        qc.assign_variants_with_valid_deletions(del_variant, {"Rv0667": locus}, {"rpoB": target})

        assert err_locus.contains_variant_with_valid_deletion(del_variant) is True
        assert err_target.contains_variant_with_valid_deletion(del_variant) is True
        assert locus.contains_variant_with_valid_deletion(del_variant) is True
        assert target.contains_variant_with_valid_deletion(del_variant) is True

    def test_deletion_outside_err_range_not_in_err(self, make_variant, make_locus_coverage, make_target_coverage):
        qc = VariantQC()
        del_variant = make_variant(nucleotide_change="c.1_100del", depth=100, freq=0.95, pos=350)

        err_locus = ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.95, average_depth=100.0)
        err_target = ERRCoverage(coords=[(310, 390)], breadth_of_coverage=0.95, average_depth=100.0)
        locus = make_locus_coverage(locus_tag="Rv0667", coords=[(100, 200), (300, 400)], err_coverage=err_locus)
        target = make_target_coverage(locus_tag="Rv0667", gene_name="rpoB", coords=[(300, 400)], err_coverage=err_target)

        qc.assign_variants_with_valid_deletions(del_variant, {"Rv0667": locus}, {"rpoB": target})

        # In locus and target (pos 350 is in range) but NOT in locus ERR (pos 350 outside ERR coords 100-200)
        assert locus.contains_variant_with_valid_deletion(del_variant) is True
        assert target.contains_variant_with_valid_deletion(del_variant) is True
        assert err_locus.contains_variant_with_valid_deletion(del_variant) is False
        # Target ERR (310-390) does contain pos 350, so the variant IS assigned to target ERR
        assert err_target.contains_variant_with_valid_deletion(del_variant) is True
