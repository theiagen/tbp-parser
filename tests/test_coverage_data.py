import pytest

class TestTargetCoverage:
    @pytest.mark.parametrize("position,expected", [
        (100, True),   # start boundary — inside
        (150, True),   # end boundary — inside
        (125, True),   # middle — inside
        (99, False),   # start - 1 — outside
        (151, False),  # end + 1 — outside
        (0, False),    # completely outside
    ])
    def test_target_coverage_contains_position(self, make_target_coverage, position, expected):
        tc = make_target_coverage(coords=(100, 150))
        assert tc.contains_position(position) == expected

    @pytest.mark.parametrize("breadth,threshold,expected", [
        (0.89, 0.90, True),   # below threshold
        (0.90, 0.90, False),  # at threshold (not strictly below)
        (0.95, 0.90, False),  # above threshold
        (0.0,  0.90, True),   # zero breadth
        (1.0,  0.90, False),  # perfect coverage
    ])
    def test_target_coverage_has_breadth_below(self, make_target_coverage, breadth, threshold, expected):
        tc = make_target_coverage(breadth_of_coverage=breadth)
        assert tc.has_breadth_below(threshold) == expected

    def test_target_coverage_breadth_above_1_raises(self, make_target_coverage):
        with pytest.raises(ValueError):
            make_target_coverage(breadth_of_coverage=1.1)


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
    def test_locus_coverage_split_contains_position(self, make_locus_coverage, position, expected):
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
    def test_locus_coverage_overlap_contains_position(self, make_locus_coverage, position, expected):
        lc = make_locus_coverage(coords=[(100, 300), (200, 400)])
        assert lc.contains_position(position) == expected

    @pytest.mark.parametrize("breadth,threshold,expected", [
        (0.89, 0.90, True),
        (0.90, 0.90, False),
        (0.95, 0.90, False),
        (0.0,  0.90, True),
        (1.0,  0.90, False),
    ])
    def test_locus_coverage_has_breadth_below(self, make_locus_coverage, breadth, threshold, expected):
        lc = make_locus_coverage(breadth_of_coverage=breadth)
        assert lc.has_breadth_below(threshold) == expected

    def test_locus_coverage_breadth_above_1_raises(self, make_locus_coverage):
        with pytest.raises(ValueError):
            make_locus_coverage(breadth_of_coverage=1.1)
