import pytest
from tbp_parser.Coverage import BedRecord, parse_bed_file

class TestBedRecordEquality:
    @pytest.mark.parametrize("field, original_value, different_value", [
            ("chrom", "chr1", "chr2"),
            ("start", 100, 400),
            ("end", 200, 500),
            ("locus_tag", "Rv0000", "Rv9999"),
            ("gene_name", "geneA",  "geneZ"),
    ])
    def test_unequal_records(self, make_bed_record, field, original_value, different_value):
        rec1 = make_bed_record(**{field: original_value})
        rec2 = make_bed_record(**{field: different_value})
        assert rec1 != rec2

    def test_records_differing_only_in_drugs_are_equal(self, make_bed_record):
        """Currently, a region is identified by its coordinates and gene, not by its drug annotations."""
        rec1 = make_bed_record(drugs=["levofloxacin"])
        rec2 = make_bed_record(drugs=["rifampicin", "isoniazid"])
        assert rec1 == rec2
        assert hash(rec1) == hash(rec2)

    def test_non_bed_record_equality(self, make_bed_record):
        rec1 = make_bed_record()
        rec2 = {"chrom": rec1.chrom, "start": rec1.start, "end": rec1.end, "locus_tag": rec1.locus_tag, "gene_name": rec1.gene_name}
        assert rec1 != rec2

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

    @pytest.mark.parametrize("alias,canonical,gene_name", [
        ("Rvnr01", "EBG00000313325", "rrs"),
        ("MTB000019", "EBG00000313325", "rrs"),
        ("Rvnr02", "EBG00000313339", "rrl"),
        ("MTB000020", "EBG00000313339", "rrl"),
        ("UNKNOWN_TAG", "UNKNOWN_TAG", "some_gene"),  # unknown kept as-is
    ], ids=["rrs-Rvnr01", "rrs-MTB000019", "rrl-Rvnr02", "rrl-MTB000020", "unknown-passthrough"])
    def test_alias_locus_tag_normalized_to_canonical(self, alias, canonical, gene_name):
        """BED locus_tag aliases should be normalized to canonical locus_tags via the GeneDatabase."""
        line = f"Chromosome\t100\t200\t{alias}\t{gene_name}"
        record = BedRecord.from_bed_line(line)
        assert record.locus_tag == canonical

    @pytest.mark.parametrize("drug_column, expected", [
        ("levofloxacin,moxifloxacin", ["levofloxacin", "moxifloxacin"]),
        ("rifampicin", ["rifampicin"]),
        (" levofloxacin , moxifloxacin ", ["levofloxacin", "moxifloxacin"]),  # whitespace stripped
        ("levofloxacin,,moxifloxacin", ["levofloxacin", "moxifloxacin"]),  # empty entries dropped
        ("fluoroquinolones,ciprofloxacin", ["fluoroquinolones", "ciprofloxacin"]),
        ("", []),
    ], ids=["two-drugs", "one-drug", "whitespace", "empty-entries", "drug-classes", "empty-column"])
    def test_drugs_parsed_from_the_sixth_column(self, drug_column, expected):
        line = f"Chromosome\t100\t200\tRv0006\tgyrA\t{drug_column}"
        record = BedRecord.from_bed_line(line)
        assert record.drugs == expected

    def test_drugs_empty_without_a_sixth_column(self):
        record = BedRecord.from_bed_line("Chromosome\t100\t200\tRv0006\tgyrA")
        assert record.drugs == []


class TestParseBedFile:
    def test_parse_records(self, tmp_path):
        bed_content = (
            "Chromosome\t100\t200\tRv0000\tgene1\n"
            "Chromosome\t300\t400\tRv0001\tgene2\n"
            "Chromosome\t500\t600\tRv0002\tgene3\n"
        )
        bed_file = tmp_path / "test.bed"
        bed_file.write_text(bed_content)

        records = parse_bed_file(str(bed_file), expected_columns=5)

        assert len(records) == 3
        assert records[0].gene_name == "gene1"
        assert records[1].locus_tag == "Rv0001"
        assert records[2].start == 500

    def test_skips_blank_lines(self, tmp_path):
        # `from_bed_line` would IndexError on [''], so blank lines have to be skipped by the parser
        bed_content = (
            "Chromosome\t100\t200\tRv0000\tgene1\n"
            "\n"
            "Chromosome\t300\t400\tRv0001\tgene2\n"
            "\n"
        )
        bed_file = tmp_path / "blank_lines.bed"
        bed_file.write_text(bed_content)

        records = parse_bed_file(str(bed_file), expected_columns=5)

        assert [record.gene_name for record in records] == ["gene1", "gene2"]

    def test_reports_every_short_line_in_a_single_error(self, tmp_path):
        # reporting them all at once means a ragged file can be fixed in one pass, and the line
        # numbers count blank lines so they line up with what an editor shows
        bed_content = (
            "Chromosome\t100\t200\tRv0000\tgene1\n"
            "Chromosome\t300\t400\n"
            "\n"
            "Chromosome\t500\t600\tRv0002\tgene3\n"
            "Chromosome\t700\n"
        )
        bed_file = tmp_path / "ragged.bed"
        bed_file.write_text(bed_content)

        with pytest.raises(ValueError, match=r"requires 5 tab-separated columns; line\(s\) 2, 5"):
            parse_bed_file(str(bed_file), expected_columns=5)

    @pytest.mark.parametrize("sixth_column", [
        "",      # no 6th column at all
        "\t",    # present but empty
        "\t ",   # present but whitespace only
    ], ids=["absent", "empty", "whitespace"])
    def test_expected_columns_of_six_requires_a_drug_column(self, tmp_path, sixth_column):
        # the `--db_bed` drug column is the truth set of gene/drug associations, so a row without one
        # cannot become a gene database entry
        bed_file = tmp_path / "db.bed"
        bed_file.write_text(f"Chromosome\t100\t200\tRv0000\tgene1{sixth_column}\n")

        with pytest.raises(ValueError, match=r"requires 6 tab-separated columns; line\(s\) 1"):
            parse_bed_file(str(bed_file), expected_columns=6)

        # the same row is a perfectly good --coverage_bed row, where the drug column is optional
        assert len(parse_bed_file(str(bed_file), expected_columns=5)) == 1


class TestUniqueBedRecordsValidation:
    def test_duplicate_records_raise_error(self, tmp_path):
        bed_content = (
            "Chromosome\t100\t200\tRv0000\tgene1\n"
            "Chromosome\t300\t400\tRv0000\tgene1\n"  # duplicate locus_tag + gene_name
        )
        bed_file = tmp_path / "test_duplicates.bed"
        bed_file.write_text(bed_content)

        with pytest.raises(ValueError, match="Duplicate BedRecords found with identical locus_tag and gene_name"):
            parse_bed_file(str(bed_file), expected_columns=5)


class TestBedRecordOverlaps:
    @pytest.mark.parametrize("b_start,b_end,expected", [
        (200, 400, True),   # left overlap
        (400, 600, False),  # no overlap (non-overlapping)
        (300, 400, True),   # adjacent: touches rec_a at position 300
        (150, 250, True),   # fully contained within rec_a
    ], ids=["left-overlap", "no-overlap", "adjacent", "full-containment"])
    def test_overlaps_with(self, make_bed_record, b_start, b_end, expected):
        a = make_bed_record(start=100, end=300, locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(start=b_start, end=b_end, locus_tag="Rv0001", gene_name="geneB")
        assert a.overlaps_with(b) is expected


class TestBedRecordOverlappingCoords:
    @pytest.mark.parametrize("a_coords,b_coords,expected", [
        ((100, 300), (200, 400), (200, 300)), # left overlap
        ((200, 400), (100, 300), (200, 300)), # right overlap
        ((100, 400), (150, 250), (150, 250)), # b fully contained in a
        ((100, 200), (200, 300), (200, 200)), # adjacent (single-point overlap)
    ], ids=["left-overlap", "right-overlap", "full-containment", "adjacent"])
    def test_overlapping_coords(self, make_bed_record, a_coords, b_coords, expected):
        a = make_bed_record(start=a_coords[0], end=a_coords[1], locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(start=b_coords[0], end=b_coords[1], locus_tag="Rv0001", gene_name="geneB")
        assert a.overlapping_coords(b) == expected

    def test_no_overlap(self, make_bed_record):
        a = make_bed_record(start=100, end=200, locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(start=300, end=400, locus_tag="Rv0001", gene_name="geneB")
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
        a = make_bed_record(start=a_coords[0], end=a_coords[1], locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(start=b_coords[0], end=b_coords[1], locus_tag="Rv0001", gene_name="geneB")
        assert a._get_non_overlapping_positions([b]) == expected_positions


class TestBedRecordGetNonOverlappingCoords:
    @pytest.mark.parametrize("a_coords,b_coords,expected_coords", [
        ((100, 200), (50, 250), []), # a fully inside b -> empty list
        ((100, 200), (150, 250), [(100, 149)]), # a=[100,200], b=[150,250] -> non-overlap = [(100, 149)]
        ((100, 400), (200, 300), [(100, 199), (301, 400)]), # a=[100,400], b=[200,300] (b inside a) -> non-overlap = [(100,199), (301,400)]
    ], ids=["fully-overlapped", "partial", "fully-encompassing"])
    def test_get_non_overlapping_coords(self, make_bed_record, a_coords, b_coords, expected_coords):
        a = make_bed_record(start=a_coords[0], end=a_coords[1], locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(start=b_coords[0], end=b_coords[1], locus_tag="Rv0001", gene_name="geneB")
        coords = a.get_non_overlapping_coords([b])
        assert sorted(coords) == sorted(expected_coords)


class TestBedRecordGetNonOverlappingReads:
    def test_reads_from_non_overlapping_positions_only(self, make_bed_record):
        # a=[100,200], b=[150,250] -> non-overlapping for a is [100,149]
        a = make_bed_record(start=100, end=200, locus_tag="Rv0000", gene_name="geneA")
        b = make_bed_record(start=150, end=250, locus_tag="Rv0001", gene_name="geneB")

        a.reads_by_position = {
            110: ["read1", "read2"],  # non-overlapping
            149: ["read3"],           # non-overlapping
            150: ["read4", "read5"],  # overlapping — excluded
            180: ["read6"],           # overlapping — excluded
        }

        reads = a.get_non_overlapping_reads([b])
        assert reads == {"read1", "read2", "read3"}
