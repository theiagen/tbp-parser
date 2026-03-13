import pytest
from tbp_parser.Coverage.coverage_data import ERRCoverage


class TestBaseCoverage:
    """Tests for methods inherited from BaseCoverage (shared across all subclasses)."""

    @pytest.mark.parametrize("breadth,threshold,expected", [
        (0.89, 0.90, True),   # below threshold
        (0.90, 0.90, False),  # at threshold (not strictly below)
        (0.95, 0.90, False),  # above threshold
        (0.0,  0.90, True),   # zero breadth
        (1.0,  0.90, False),  # perfect coverage
    ])
    def test_has_breadth_below(self, make_locus_coverage, breadth, threshold, expected):
        lc = make_locus_coverage(breadth_of_coverage=breadth)
        assert lc.has_breadth_below(threshold) == expected

    def test_breadth_above_1_raises(self, make_locus_coverage):
        with pytest.raises(ValueError):
            make_locus_coverage(breadth_of_coverage=1.1)


class TestTargetCoverage:
    @pytest.mark.parametrize("position,expected", [
        (100, True),   # start boundary — inside
        (150, True),   # end boundary — inside
        (125, True),   # middle — inside
        (99, False),   # start - 1 — outside
        (151, False),  # end + 1 — outside
        (0, False),    # completely outside
    ])
    def test_contains_position(self, make_target_coverage, position, expected):
        tc = make_target_coverage(coords=[(100, 150)])
        assert tc.contains_position(position) == expected


class TestLocusCoverage:
    @pytest.mark.parametrize("position,expected", [
        (100, True),   # start of range 1
        (150, True),   # end of range 1
        (125, True),   # middle of range 1
        (250, True),   # start of range 2
        (350, True),   # end of range 2
        (300, True),   # middle of range 2
        (200, False),  # between the two ranges
        (50,  False),  # before range 1
        (400, False),  # after range 2
    ])
    def test_split_ranges_contains_position(self, make_locus_coverage, position, expected):
        lc = make_locus_coverage(coords=[(100, 150), (250, 350)])
        assert lc.contains_position(position) == expected

    @pytest.mark.parametrize("position,expected", [
        (100, True),   # start of range 1
        (150, True),   # middle of range 1
        (300, True),   # end of range 1
        (200, True),   # start of range 2
        (350, True),   # middle of range 2
        (400, True),   # end of range 2
        (50,  False),  # before range 1
        (450, False),  # after range 2
    ])
    def test_overlapping_ranges_contains_position(self, make_locus_coverage, position, expected):
        lc = make_locus_coverage(coords=[(100, 300), (200, 400)])
        assert lc.contains_position(position) == expected


class TestOverlapsRange:
    """Tests for overlaps_range on TargetCoverage and LocusCoverage."""

    @pytest.mark.parametrize("start,end,expected", [
        (100, 150, True),   # exact match
        (90, 110, True),    # partial overlap at start
        (140, 160, True),   # partial overlap at end
        (80, 200, True),    # fully encompassing
        (110, 140, True),   # fully contained
        (151, 200, False),  # adjacent but outside (just past end)
        (50, 99, False),    # adjacent but outside (just before start)
        (200, 300, False),  # completely outside
    ])
    def test_target_overlaps_range(self, make_target_coverage, start, end, expected):
        tc = make_target_coverage(coords=[(100, 150)])
        assert tc.overlaps_range(start, end) == expected

    @pytest.mark.parametrize("start,end,expected", [
        (100, 150, True),   # exact match range 1
        (250, 350, True),   # exact match range 2
        (90, 110, True),    # partial overlap at start of range 1
        (340, 360, True),   # partial overlap at end of range 2
        (80, 400, True),    # fully encompassing both ranges
        (120, 130, True),   # fully contained in range 1
        (151, 249, False),  # between the two ranges
        (50, 99, False),    # before range 1
        (351, 400, False),  # after range 2
    ])
    def test_locus_overlaps_range(self, make_locus_coverage, start, end, expected):
        lc = make_locus_coverage(coords=[(100, 150), (250, 350)])
        assert lc.overlaps_range(start, end) == expected

    def test_variant_deletion_spanning_boundary(self, make_variant, make_target_coverage):
        """A deletion starting outside but spanning into the coverage region uses overlaps_range."""
        # c.1_20del with pos=90 -> absolute_start=90, absolute_end=90+19=109
        v = make_variant(nucleotide_change="c.1_20del", pos=90)
        tc = make_target_coverage(coords=[(100, 150)])
        # Point check fails (90 < 100), but range overlap succeeds
        assert tc.contains_position(v.pos) is False
        assert tc.overlaps_range(v.absolute_start, v.absolute_end) is True

    def test_variant_point_mutation_inside(self, make_variant, make_target_coverage):
        """A SNP inside the region (absolute_start == absolute_end == pos)."""
        v = make_variant(nucleotide_change="c.1349C>T", pos=120)
        tc = make_target_coverage(coords=[(100, 150)])
        assert tc.overlaps_range(v.absolute_start, v.absolute_end) is True

    def test_variant_deletion_fully_outside(self, make_variant, make_target_coverage):
        """A deletion that doesn't reach the coverage region."""
        # c.1_5del with pos=50 -> absolute_start=50, absolute_end=54
        v = make_variant(nucleotide_change="c.1_5del", pos=50)
        tc = make_target_coverage(coords=[(100, 150)])
        assert tc.overlaps_range(v.absolute_start, v.absolute_end) is False


class TestERRWithinCoords:
    """Test that Locus/TargetCoverage with an attached ERRCoverage validates that the ERR coords are within bounds."""
    @pytest.mark.parametrize("err_coords", [
        [(140, 160)], # fully within target/locus
        [(100, 200)], # exactly matches target/locus (not strictly within)
    ])
    def test_err_within_coords_valid(self, make_target_coverage, make_locus_coverage, err_coords):
        err = ERRCoverage(coords=err_coords, breadth_of_coverage=0.95, average_depth=50.0)
        tc = make_target_coverage(coords=[(100, 200)], err_coverage=err)
        lc = make_locus_coverage(coords=[(100, 200)], err_coverage=err)
        assert tc.err_coverage == err
        assert lc.err_coverage == err

    @pytest.mark.parametrize("err_coords", [
        [(90, 110)], # starts before target/locus
        [(190, 210)], # ends after target/locus
        [(200, 210)], # starts exactly at end boundary but extends beyond
        [(90, 100)],  # ends exactly at start boundary but starts before
        [(50, 250)] # completely encompasses target/locus
    ])
    def test_err_within_coords_invalid(self, make_target_coverage, make_locus_coverage, err_coords):
        err = ERRCoverage(coords=err_coords, breadth_of_coverage=0.95, average_depth=50.0)
        with pytest.raises(ValueError):
            make_target_coverage(coords=[(100, 200)], err_coverage=err)

        with pytest.raises(ValueError):
            make_locus_coverage(coords=[(100, 200)], err_coverage=err)
class TestContainsValidDeletion:
    def test_contains_valid_deletion(self, make_variant):
        del_variant = make_variant(nucleotide_change="c.1_100del", pos=150)
        err = ERRCoverage(coords=[(100, 200)], breadth_of_coverage=0.95, average_depth=50.0, valid_deletions=[del_variant])
        assert err.contains_valid_deletion(del_variant) is True
