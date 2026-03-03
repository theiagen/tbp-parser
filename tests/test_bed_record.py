import pytest
from Coverage import BedRecord, parse_bed_file

class TestBedRecordFromBedLine:
    def test_basic_parsing(self):
        line = "Chromosome\t759807\t763325\tRv0667\trpoB"
        record = BedRecord.from_bed_line(line)
        assert record.chrom == "Chromosome"
        assert record.start == 759807
        assert record.end == 763325
        assert record.locus_tag == "Rv0667"
        assert record.gene_name == "rpoB"

    def test_derived_fields(self):
        line = "Chromosome\t100\t200\tRv0000\tgene1"
        record = BedRecord.from_bed_line(line)
        assert record.length == 101        # end - start + 1 (1-based)
        assert record.coords == (100, 200)

    @pytest.mark.parametrize("input_gene_name, expected_gene_name", [
        ("mmpR5", "Rv0678"),
        ("fbiD", "Rv2983"),
    ])
    def test_gene_name_normalization(self, input_gene_name, expected_gene_name):
        """mmpR5 is normalized to Rv0678 by Helper.normalize_field_values."""
        line = f"Chromosome\t100\t200\tRv0678\t{input_gene_name}"
        record = BedRecord.from_bed_line(line)
        assert record.gene_name == expected_gene_name


class TestParseBedFile:
    def test_parse_records(self, tmp_path):
        bed_content = (
            "Chromosome\t100\t200\tRv0000\tgene1\n"
            "Chromosome\t300\t400\tRv0001\tgene2\n"
            "Chromosome\t500\t600\tRv0002\tgene3\n"
        )
        bed_file = tmp_path / "test.bed"
        bed_file.write_text(bed_content)

        records = parse_bed_file(str(bed_file))

        assert len(records) == 3
        assert records[0].gene_name == "gene1"
        assert records[1].locus_tag == "Rv0001"
        assert records[2].start == 500


class TestUniqueBedRecordsValidation:
    def test_duplicate_records_raise_error(self, tmp_path):
        bed_content = (
            "Chromosome\t100\t200\tRv0000\tgene1\n"
            "Chromosome\t300\t400\tRv0000\tgene1\n"  # duplicate locus_tag + gene_name
        )
        bed_file = tmp_path / "test_duplicates.bed"
        bed_file.write_text(bed_content)

        with pytest.raises(ValueError, match="Duplicate BedRecords found with identical locus_tag and gene_name"):
            parse_bed_file(str(bed_file))


class TestBedRecordOverlaps:
    @pytest.mark.parametrize("b_start,b_end,expected", [
        (200, 400, True),   # left overlap
        (400, 600, False),  # no overlap (non-overlapping)
        (300, 400, True),   # adjacent: touches rec_a at position 300
        (150, 250, True),   # fully contained within rec_a
    ], ids=["left-overlap", "no-overlap", "adjacent", "full-containment"])
    def test_overlaps_with(self, make_bed_record, b_start, b_end, expected):
        a = make_bed_record(100, 300, locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(b_start, b_end, locus_tag="Rv0001", gene_name="geneB")
        assert a.overlaps_with(b) is expected


class TestBedRecordOverlappingCoords:
    @pytest.mark.parametrize("a_coords,b_coords,expected", [
        ((100, 300), (200, 400), (200, 300)), # left overlap
        ((200, 400), (100, 300), (200, 300)), # right overlap
        ((100, 400), (150, 250), (150, 250)), # b fully contained in a
        ((100, 200), (200, 300), (200, 200)), # adjacent (single-point overlap)
    ], ids=["left-overlap", "right-overlap", "full-containment", "adjacent"])
    def test_overlapping_coords(self, make_bed_record, a_coords, b_coords, expected):
        a = make_bed_record(a_coords[0], a_coords[1], locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(b_coords[0], b_coords[1], locus_tag="Rv0001", gene_name="geneB")
        assert a.overlapping_coords(b) == expected

    def test_no_overlap(self, make_bed_record):
        a = make_bed_record(100, 200, locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(300, 400, locus_tag="Rv0001", gene_name="geneB")
        with pytest.raises(Exception, match="No overlap"):
            a.overlapping_coords(b)


class TestBedRecordGetNonOverlappingPositions:
    @pytest.mark.parametrize("a_coords,b_coords,expected_positions", [
        ((100, 200), (150, 250), set(range(100, 150))), # a=[100,200], b=[150,250] -> non-overlapping for a = [100,149]
        ((100, 200), (50, 150), set(range(151, 201))), # a=[100,200], b=[50,150] -> non-overlapping for a = [151,200]
        ((100, 200), (50, 250), set()), # a fully inside b -> no non-overlapping positions
        ((50, 250), (100, 200), set(range(50, 100)) | set(range(201, 251))), # a fully encompasses b -> non-overlapping = [50,99] ∪ [201,250]
    ], ids=["left-overlap", "right-overlap", "full-containment", "fully-encompassing"])
    def test_get_non_overlapping_positions(self, make_bed_record, a_coords, b_coords, expected_positions):
        a = make_bed_record(a_coords[0], a_coords[1], locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(b_coords[0], b_coords[1], locus_tag="Rv0001", gene_name="geneB")
        assert a._get_non_overlapping_positions([b]) == expected_positions


class TestBedRecordGetNonOverlappingCoords:
    @pytest.mark.parametrize("a_coords,b_coords,expected_coords", [
        ((100, 200), (50, 250), []), # a fully inside b -> empty list
        ((100, 200), (150, 250), [(100, 149)]), # a=[100,200], b=[150,250] -> non-overlap = [(100, 149)]
        ((100, 400), (200, 300), [(100, 199), (301, 400)]), # a=[100,400], b=[200,300] (b inside a) -> non-overlap = [(100,199), (301,400)]
    ], ids=["fully-overlapped", "partial", "fully-encompassing"])
    def test_get_non_overlapping_coords(self, make_bed_record, a_coords, b_coords, expected_coords):
        a = make_bed_record(a_coords[0], a_coords[1], locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(b_coords[0], b_coords[1], locus_tag="Rv0001", gene_name="geneB")
        coords = a.get_non_overlapping_coords([b])
        assert sorted(coords) == sorted(expected_coords)


class TestBedRecordGetNonOverlappingReads:
    def test_reads_from_non_overlapping_positions_only(self, make_bed_record):
        # a=[100,200], b=[150,250] -> non-overlapping for a is [100,149]
        a = make_bed_record(100, 200, locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(150, 250, locus_tag="Rv0001", gene_name="geneB")

        a.reads_by_position = {
            110: ["read1", "read2"],  # non-overlapping
            149: ["read3"],           # non-overlapping
            150: ["read4", "read5"],  # overlapping — excluded
            180: ["read6"],           # overlapping — excluded
        }

        reads = a.get_non_overlapping_reads([b])
        assert reads == {"read1", "read2", "read3"}
