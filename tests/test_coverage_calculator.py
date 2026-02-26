from collections import defaultdict

class TestPopulateReadsByPosition:
    def test_populate_reads_by_position(self, make_cov_calc, make_bed_record):
        calc = make_cov_calc(cov_start=0, cov_end=100, read_length=10)
        self.config = calc.config

        record1 = make_bed_record(1, 50, locus_tag="Rv0000", gene_name="geneA")
        calc.populate_reads_by_position([record1])
        expected = {
            1: ["read1", "read1"], # counting mates
            2: ["read1", "read1", "read2", "read2"],
            10: [
                "read1", "read1", "read2", "read2",
                "read3", "read3", "read4", "read4",
                "read5", "read5", "read6", "read6",
                "read7", "read7", "read8", "read8",
                "read9", "read9", "read10", "read10",
            ],
            50: [
                "read41", "read41", "read42", "read42",
                "read43", "read43", "read44", "read44",
                "read45", "read45", "read46", "read46",
                "read47", "read47", "read48", "read48",
                "read49", "read49", "read50", "read50",
            ]
        }
        for pos, reads in expected.items():
            assert record1.reads_by_position[pos] == reads


class TestResolveOverlappingRegions:
    def test_resolve_overlapping_regions(self, make_cov_calc, make_bed_record):
        calc = make_cov_calc(cov_start=0, cov_end=100, read_length=10)
        self.config = calc.config

        record1 = make_bed_record(10, 60, locus_tag="Rv0000", gene_name="geneA")
        record2 = make_bed_record(40, 90, locus_tag="Rv0000", gene_name="geneB")
        all_records = [record1, record2]

        calc.populate_reads_by_position(all_records)

        # before resolving overlaps confirm that both records have reads in the overlapping region
        assert record1.overlaps_with(record2) == True
        assert record2.overlaps_with(record1) == True

        assert record1.get_non_overlapping_coords([record2]) == [(10, 39)]
        assert record2.get_non_overlapping_coords([record1]) == [(61, 90)]
        # getting non-overlapping reads compared to nothing returns all unique reads for that record
        # proving overlapping reads exist in both records
        assert record1.get_non_overlapping_reads([]) & set(f"read{i}" for i in range(40, 61))
        assert record2.get_non_overlapping_reads([]) & set(f"read{i}" for i in range(40, 61))

        # after resolving overlaps, the reads in the overlapping region should be de-duplicated across the two records
        # reads existing in the non-overlapping region can extend into the overlapping region which is why we only don't see reads 40-51 (11 total)
        calc.resolve_overlapping_regions(all_records)
        assert not record1.get_non_overlapping_reads([record2]) & set(f"read{i}" for i in range(40, 52))
        assert not record2.get_non_overlapping_reads([record1]) & set(f"read{i}" for i in range(40, 52))


class TestCalculateBreadthOfCoverage:
    def test_all_above_min_depth(self, make_cov_calc):
        calc = make_cov_calc()
        calc.config.MIN_DEPTH = 10

        reads_by_position = {i: [f"r{j}" for j in range(10)] for i in range(100, 110)}
        result = calc._calculate_breadth_of_coverage(reads_by_position)
        assert result == 1.0

    def test_some_below_min_depth(self, make_cov_calc):
        calc = make_cov_calc()
        calc.config.MIN_DEPTH = 15

        reads_by_position = {
            1: ["r"] * 15,  # at min depth -> counts
            2: ["r"] * 16,  # above -> counts
            3: ["r"] * 5,   # below -> does not count
            4: ["r"] * 2,   # below -> does not count
            5: ["r"] * 0,  # below -> does not count
        }
        result = calc._calculate_breadth_of_coverage(reads_by_position)
        assert result == 0.4  # 2 out of 5

    def test_none_above_min_depth(self, make_cov_calc):
        calc = make_cov_calc()
        reads_by_position = {1: [], 2: ["r"] * 2, 3: []}
        result = calc._calculate_breadth_of_coverage(reads_by_position)
        assert result == 0.0

    def test_zero_positions(self, make_cov_calc):
        calc = make_cov_calc()
        result = calc._calculate_breadth_of_coverage({})
        assert result == 0.0


class TestCalculateAverageDepth:
    def test_basic_average(self, make_cov_calc):
        calc = make_cov_calc()
        reads_by_position = {
            1: ["r1", "r2", "r3"],  # depth 3
            2: ["r1", "r2"],        # depth 2
            3: ["r1"],              # depth 1
        }
        result = calc._calculate_average_depth(reads_by_position)
        assert result == 2.0  # (3+2+1) / 3

    def test_zero_positions(self, make_cov_calc):
        calc = make_cov_calc()
        result = calc._calculate_average_depth({})
        assert result == 0.0


class TestGenerateCoverageMaps:
    def test_empty_reads_gives_zero_breadth(self, make_cov_calc, make_bed_record):
        calc = make_cov_calc(cov_start=50, cov_end=100, read_length=10)
        calc.config.MIN_DEPTH = 0
        self.config = calc.config

        # coverage only covers 50-100, so reads_by_position will be empty.
        record = make_bed_record(5, 20, locus_tag="Rv0000", gene_name="geneA")
        target_map, locus_map = calc.generate_coverage_maps([record])
        assert target_map["geneA"].breadth_of_coverage == 0.0
        assert target_map["geneA"].average_depth == 0.0
        assert locus_map["Rv0000"].breadth_of_coverage == 0.0
        assert locus_map["Rv0000"].average_depth == 0.0

    def test_single_record_below_min_depth(self, make_cov_calc, make_bed_record):
        calc = make_cov_calc(cov_start=0, cov_end=100, read_length=10)
        calc.config.MIN_DEPTH = 100
        self.config = calc.config

        # each position from 50-50 should have 10 read pairs,
        record = make_bed_record(50, 60, locus_tag="Rv0000", gene_name="geneA")
        calc.populate_reads_by_position([record])

        target_map, locus_map = calc.generate_coverage_maps([record])
        assert target_map["geneA"].breadth_of_coverage == 0.0
        assert target_map["geneA"].average_depth == 20.0
        assert locus_map["Rv0000"].breadth_of_coverage == 0.0
        assert locus_map["Rv0000"].average_depth == 20.0

    def test_single_record_above_min_depth(self, make_cov_calc, make_bed_record):
        calc = make_cov_calc(cov_start=0, cov_end=100, read_length=10)
        calc.config.MIN_DEPTH = 5
        self.config = calc.config

        # each position from 50-50 should have 10 read pairs,
        record = make_bed_record(50, 60, locus_tag="Rv0000", gene_name="geneA")
        calc.populate_reads_by_position([record])

        target_map, locus_map = calc.generate_coverage_maps([record])
        assert target_map["geneA"].breadth_of_coverage == 1.0
        assert target_map["geneA"].average_depth == 20.0
        assert locus_map["Rv0000"].breadth_of_coverage == 1.0
        assert locus_map["Rv0000"].average_depth == 20.0

class TestComplexOverlappingRecords:
    def test_complex_overlapping_records(self, make_cov_calc, make_bed_record):
        calc = make_cov_calc(cov_start=0, cov_end=100, read_length=10)
        calc.config.MIN_DEPTH = 8
        self.config = calc.config

        record1 = make_bed_record(1, 40, locus_tag="Rv0000", gene_name="geneA")
        record2 = make_bed_record(25, 75, locus_tag="Rv0000", gene_name="geneB")
        record3 = make_bed_record(55, 100, locus_tag="Rv0000", gene_name="geneC")
        all_records = [record1, record2, record3]

        calc.populate_reads_by_position(all_records)

        # before resolving overlaps confirm that all records overlap (if applicable)
        assert record1.overlaps_with(record2) == True
        assert record2.overlaps_with(record1) == True
        assert record2.overlaps_with(record3) == True
        assert record3.overlaps_with(record2) == True
        assert record1.overlaps_with(record3) == False

        # before resolving overlaps, the non-overlapping coords should be as follows:
        assert record1.get_non_overlapping_coords([record2, record3]) == [(1, 24)]
        assert record2.get_non_overlapping_coords([record1, record3]) == [(41, 54)]
        assert record3.get_non_overlapping_coords([record1, record2]) == [(76, 100)]

        # Coverage increases linearly from 1 to 10 read-pairs across positions 1–10,
        # giving an average depth of 11 and 0.7 breadth (positions 1–3 < 8 reads).
        first_10_positions = {k: v for k, v in record1.reads_by_position.items() if k <= 10}
        assert calc._calculate_breadth_of_coverage(first_10_positions) == 0.7
        assert calc._calculate_average_depth(first_10_positions) == 11.0

        # before resolving overlaps, the reads in the overlapping regions should be as follows:
        overlap_coords = [
            record1.overlapping_coords(record2),
            record2.overlapping_coords(record1),
            record2.overlapping_coords(record3),
            record3.overlapping_coords(record2),
        ]
        assert overlap_coords == [(25, 40), (25, 40), (55, 75), (55, 75)]

        # confirming that overlapping regions have inflated depths
        record1_overlap_positions = {
            k: v for k, v in record1.reads_by_position.items()
            for k in range(25, 41)
        }
        record2_overlap_positions = {
            k: v for k, v in record2.reads_by_position.items()
            for k in list(range(25, 41)) + list(range(55, 76))
        }
        record3_overlap_positions = {
            k: v for k, v in record3.reads_by_position.items()
            for k in range(55, 76)
        }
        all_overlap_reads_by_position = defaultdict(list)
        for d in [record1_overlap_positions, record2_overlap_positions, record3_overlap_positions]:
            for pos, read_list in d.items():
                all_overlap_reads_by_position[pos].extend(read_list)

        # each position would normally have 10 read pairs, but this is doubled in the overlapping regions (40 total)
        assert calc._calculate_breadth_of_coverage(all_overlap_reads_by_position) == 1.0
        assert calc._calculate_average_depth(all_overlap_reads_by_position) == 40.0

        # before resolving overlaps, the breadth of coverage and average depth should be as follows:
        target_map, locus_map = calc.generate_coverage_maps(all_records)

        # BedRecord([geneA][Rv0000](1, 40))
        assert target_map["geneA"].breadth_of_coverage == 0.925 # 37/40 >= min depth
        assert target_map["geneA"].average_depth == 17.75 # 710/40

        # BedRecord([geneB][Rv0000](25, 75))
        assert target_map["geneB"].breadth_of_coverage == 1.0 # 51/51 >= min depth
        assert target_map["geneB"].average_depth == 20.0 # 1020/51

        # BedRecord([geneC][Rv0000](55, 100))
        assert target_map["geneC"].breadth_of_coverage == 1.0 # 46/46 >= min depth
        assert target_map["geneC"].average_depth == 20.0 # 920/46

        # AVG DEPTH ACROSS LOCUS:
        # [1-10]: 11 (110 reads) (10 pos)
        # [11-24]: 20 (280 reads) (14 pos)
        # [25-40]: 20*2 (320*2 reads) (16 pos)
        # [41-54]: 20 (280 reads) (14 pos)
        # [55-75]: 20*2 (420*2 reads) (21 pos)
        # [76-100]: 20 (500 reads) (25 pos)
        assert locus_map["Rv0000"].breadth_of_coverage == 0.97 # 97/100 >= min depth
        assert locus_map["Rv0000"].average_depth == 26.50 # 2650/100