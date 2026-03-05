import pytest
from Variant import Variant, VariantRecord, VariantProcessor
from Utilities import GeneDatabase

class TestExpandConsequences:
    """Tests for VariantProcessor._expand_consequences"""

    def test_invalid_consequences_returns_only_original(self, make_variant_record):
        """Non [mmpR5, mmpL5, mmpS5] gene should return only the original record."""
        vr = make_variant_record(gene_id="Rv0000", gene_name="geneA")  # Not mmpR5/mmpL5/mmpS5
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 1
        assert result[0] is vr

    def test_valid_empty_consequences_returns_original(self, make_variant_record):
        """Rv0678 with empty consequences should return only original."""
        vr = make_variant_record(gene_id="Rv0678", consequences=[])
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 1
        assert result[0] is vr

    def test_valid_consequences_Rv0676c_expansion(self, make_variant_record, make_consequences):
        """Rv0676c with valid consequences should create 2 records."""
        valid_csq = make_consequences(gene_id="Rv0678")
        vr = make_variant_record(gene_id="Rv0676c", consequences=[valid_csq])
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 2
        assert result[0] is vr
        assert result[0].gene_id == "Rv0676c"
        assert type(result[1]) is VariantRecord
        assert result[1].gene_id == "Rv0678"

    def test_valid_consequences_Rv0677c_expansion(self, make_variant_record, make_consequences):
        """Rv0677c with valid consequences should create 2 records."""
        valid_csq = make_consequences(gene_id="Rv0678")
        vr = make_variant_record(gene_id="Rv0677c", consequences=[valid_csq])
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 2
        assert result[0] is vr
        assert result[0].gene_id == "Rv0677c"
        assert type(result[1]) is VariantRecord
        assert result[1].gene_id == "Rv0678"

    def test_valid_consequences_Rv0678_expansion(self, make_variant_record, make_consequences):
        """Rv0678 with valid consequences should create 2 records."""
        valid_csq = make_consequences(gene_id="Rv0676c")
        vr = make_variant_record(gene_id="Rv0678", consequences=[valid_csq])
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 2
        assert result[0] is vr
        assert result[0].gene_id == "Rv0678"
        assert type(result[1]) is VariantRecord
        assert result[1].gene_id == "Rv0676c"

    def test_valid_consequences_multiple_expansion(self, make_variant_record, make_consequences):
        """Rv0676c with valid consequences should create 3 records.
        Consequences with same gene_id as parent should be skipped."""
        mmpL5_csq = make_consequences(gene_id="Rv0676c", gene_name="skip_same_parent")
        mmpS5_csq = make_consequences(gene_id="Rv0677c", gene_name="mmpS5")
        mmpR5_csq = make_consequences(gene_id="Rv0678", gene_name="mmpR5")

        vr = make_variant_record(
            gene_id="Rv0676c",
            gene_name="mmpL5",
            consequences=[
                mmpL5_csq,
                mmpS5_csq,
                mmpR5_csq,
            ]
        )
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 3
        gene_names = [r.gene_name for r in result]
        # note gene_name != "skip_same_parent"
        assert gene_names == ["mmpL5", "mmpS5", "mmpR5"]

    def test_valid_consequences_multiple_expansion_random_genes(self, make_variant_record, make_consequences):
        """
        NOTE: Not sure if this should be intended behavior or not
        Rv0676c with valid consequences containing mix of genes should still create multiple records.
        Consequences with same gene_id as parent should still be skipped.
        """
        mmpL5_csq = make_consequences(gene_id="Rv0676c", gene_name="skip_same_parent")
        mmpS5_csq = make_consequences(gene_id="Rv0677c", gene_name="mmpS5")
        mmpR5_csq = make_consequences(gene_id="Rv0678", gene_name="mmpR5")
        random_csq = make_consequences(gene_id="Rv0000", gene_name="geneA")

        vr = make_variant_record(
            gene_id="Rv0676c",
            gene_name="mmpL5",
            consequences=[
                mmpL5_csq,
                mmpS5_csq,
                mmpR5_csq,
                random_csq,
            ]
        )
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 4
        gene_names = [r.gene_name for r in result]
        # note gene_name != "skip_same_parent"
        assert gene_names == ["mmpL5", "mmpS5", "mmpR5", "geneA"]

    def test_valid_consequences_overwrite_vr_fields(self, make_variant_record, make_annotation, make_consequences):
        """Expanded record should just preserve pos/depth/freq from parent and inherit everything else from consequence/annotation."""
        consequence = make_consequences(
            gene_id="Rv0676c",
            gene_name="mmpL5",
            type="upstream_gene_variant",
            nucleotide_change="c.-200C>T",
            protein_change="p.Val830Ile",
            annotation=[
                make_annotation(
                    drug="clofazimine",
                    confidence="Uncertain significance",
                    source="foo",
                    comment="bar",
                )
            ],
        )

        vr = make_variant_record(
            gene_id="Rv0678",
            gene_name="mmpR5",
            pos=778990,
            depth=200,
            freq=0.80,
            nucleotide_change="c.192C>T",
            protein_change="p.Asp64Asn",
            annotation=[
                make_annotation(
                    drug="bedaquiline",
                    confidence="Assoc w R",
                    source="",
                    comment="",
                )
            ],
            consequences=[consequence],
        )
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        csq_vr = result[1]
        assert csq_vr.gene_id == "Rv0676c"
        assert csq_vr.gene_name == "mmpL5"
        assert csq_vr.pos == 778990
        assert csq_vr.depth == 200
        assert csq_vr.freq == 0.80
        assert csq_vr.nucleotide_change == "c.-200C>T"
        assert csq_vr.protein_change == "p.Val830Ile"
        assert len(csq_vr.annotation) == 1
        assert csq_vr.annotation[0].drug == "clofazimine"
        assert csq_vr.annotation[0].confidence == "Uncertain significance"
        assert csq_vr.annotation[0].source == "foo"
        assert csq_vr.annotation[0].comment == "bar"

    def test_valid_consequences_overwrite_vr_fields_no_annotation(self, make_variant_record, make_annotation, make_consequences):
        """Expanded record should just preserve pos/depth/freq from parent and inherit everything else from consequence/annotation.
        If consequence has no annotation, it still expands, but parses gene_associated_drugs field for that information.
        See variant_record.from_consequences method."""
        consequence = make_consequences(
            gene_id="Rv0676c",
            gene_name="mmpL5",
            type="upstream_gene_variant",
            nucleotide_change="c.-200C>T",
            protein_change="p.Val830Ile",
            annotation=[],
        )

        vr = make_variant_record(
            gene_id="Rv0678",
            gene_name="mmpR5",
            pos=778990,
            depth=200,
            freq=0.80,
            nucleotide_change="c.192C>T",
            protein_change="p.Asp64Asn",
            annotation=[
                make_annotation(
                    drug="bedaquiline",
                    confidence="Assoc w R",
                    source="foo",
                    comment="bar",
                )
            ],
            consequences=[consequence],
            gene_associated_drugs=["drug1", "drug2"],
        )
        processor = VariantProcessor()
        result = processor._expand_consequences(vr)
        assert len(result) == 2

        csq_vr = result[1]
        assert csq_vr.gene_id == "Rv0676c"
        assert csq_vr.gene_name == "mmpL5"
        assert csq_vr.pos == 778990
        assert csq_vr.depth == 200
        assert csq_vr.freq == 0.80
        assert csq_vr.nucleotide_change == "c.-200C>T"
        assert csq_vr.protein_change == "p.Val830Ile"
        assert len(csq_vr.annotation) == 2
        assert csq_vr.annotation[0].drug == "drug1"
        assert csq_vr.annotation[0].confidence == "No WHO annotation"
        assert csq_vr.annotation[0].source == ""
        assert csq_vr.annotation[0].comment == ""
        assert csq_vr.annotation[1].drug == "drug2"
        assert csq_vr.annotation[1].confidence == "No WHO annotation"
        assert csq_vr.annotation[1].source == ""
        assert csq_vr.annotation[1].comment == ""


class TestGetVariantsFromAnnotations:
    """Tests for VariantProcessor._get_variants_from_annotations"""

    def test_empty_annotations_returns_empty(self, make_variant_record):
        """Record with empty annotations should return empty list."""
        vr = make_variant_record(annotation=[])
        processor = VariantProcessor()
        result = processor._get_variants_from_annotations(vr)
        assert result == []

    def test_single_annotation_creates_variant(self, make_variant_record, make_annotation):
        """
        Single annotation should create Variants based on:
        - drugs in the annotation field
        - drugs in the gene_associated_drugs field
        - drugs in the GeneDatabase for that gene
        """
        anno = make_annotation(drug="foo", confidence="Assoc w R")

        # drugs associated with Rv0005 (gyrB) in GeneDatabase are: ["levofloxacin", "moxifloxacin"]
        vr = make_variant_record(
            gene_id="Rv0005",
            gene_name="gyrB",
            annotation=[anno],
            gene_associated_drugs=["rifampicin"],
        )
        processor = VariantProcessor()
        result = processor._get_variants_from_annotations(vr)
        result = sorted(result, key=lambda v: v.drug)
        assert len(result) == 4
        assert [type(v) for v in result] == [Variant, Variant, Variant, Variant]
        assert [v.drug for v in result] == ["foo", "levofloxacin", "moxifloxacin", "rifampicin"]

    def test_gene_database_drugs_create_synthetic_variants(self, make_variant_record, make_annotation):
        """Drugs from GeneDatabase not yet seen should create synthetic variants."""
        # gyrB (Rv0005): GeneDatabase drugs = ["levofloxacin", "moxifloxacin"]
        anno = make_annotation(drug="levofloxacin", confidence="Assoc w R")
        vr = make_variant_record(
            gene_id="Rv0005",
            gene_name="gyrB",
            nucleotide_change="c.1500A>G",
            protein_change="p.Ala500Thr",
            annotation=[anno],
            gene_associated_drugs=["levofloxacin"],
        )
        processor = VariantProcessor()
        result = processor._get_variants_from_annotations(vr)
        drugs = {v.drug for v in result}
        assert "moxifloxacin" in drugs

    def test_multiple_annotations_create_multiple_variants(self, make_variant_record, make_annotation):
        """Multiple annotations should create multiple variants."""
        anno1 = make_annotation(drug="foo", confidence="Assoc w R")
        anno2 = make_annotation(drug="bar", confidence="Assoc w R")
        vr = make_variant_record(
            gene_id="Rv0678",
            gene_name="mmpR5",
            nucleotide_change="c.192C>T",
            protein_change="p.Asp64Asn",
            annotation=[anno1, anno2],
            gene_associated_drugs=["bedaquiline", "clofazimine"],
        )
        processor = VariantProcessor()
        result = processor._get_variants_from_annotations(vr)
        drugs = {v.drug for v in result}
        assert "foo" in drugs
        assert "bar" in drugs


class TestExpandAnnotationsForAllDrugs:
    """Tests for VariantProcessor._expand_annotations_for_all_drugs"""

    def test_original_annotations_in_result(self, make_annotation):
        """Original annotations should be in the result set."""
        anno = make_annotation(drug="rifampin", confidence="Assoc w R", source="WHO")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=["rifampin"],
            gene_id="Rv0667",
        )
        assert anno in result

    def test_gene_associated_drugs_added(self, make_annotation):
        """Gene-associated drugs not in original should be added with 'No WHO annotation'."""
        anno = make_annotation(drug="bedaquiline", confidence="Assoc w R")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=["bedaquiline", "clofazimine"],
            gene_id="Rv0678",
        )
        drugs = {a.drug for a in result}
        assert "clofazimine" in drugs
        clofazimine_annos = [a for a in result if a.drug == "clofazimine"]
        assert all(a.confidence == "No WHO annotation" for a in clofazimine_annos)

    def test_gene_database_drugs_added(self, make_annotation):
        """GeneDatabase drugs not already seen should be added with 'not in TBDB' source."""
        # gyrB (Rv0005): GeneDatabase has ["levofloxacin", "moxifloxacin"]
        anno = make_annotation(drug="levofloxacin", confidence="Assoc w R")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=["levofloxacin"],
            gene_id="Rv0005",
        )
        drugs = {a.drug for a in result}
        assert "moxifloxacin" in drugs
        moxi_annos = [a for a in result if a.drug == "moxifloxacin"]
        assert any("not in TBDB" in a.source for a in moxi_annos)

    def test_empty_gene_associated_drugs_still_adds_gene_database(self, make_annotation):
        """Empty gene_associated_drugs should still add GeneDatabase drugs."""
        anno = make_annotation(drug="levofloxacin", confidence="Assoc w R")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=[],
            gene_id="Rv0005",  # GeneDatabase: ["levofloxacin", "moxifloxacin"]
        )
        drugs = {a.drug for a in result}
        assert "moxifloxacin" in drugs

    def test_gene_associated_drug_already_in_original_not_duplicated(self, make_annotation):
        """Gene-associated drug already in original should not create extra annotation."""
        anno = make_annotation(drug="rifampin", confidence="Assoc w R", source="WHO")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=["rifampin"],
            gene_id="Rv0667",
        )
        rifampin_annos = [a for a in result if a.drug == "rifampin"]
        assert len(rifampin_annos) == 1
        assert rifampin_annos[0].source == "WHO"

    def test_gene_database_drug_already_in_gene_associated_not_duplicated(self, make_annotation):
        """GeneDatabase drug already covered by gene_associated_drugs should not be re-added."""
        anno = make_annotation(drug="levofloxacin", confidence="Assoc w R")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=["levofloxacin", "moxifloxacin"],
            gene_id="Rv0005",  # GeneDatabase: ["levofloxacin", "moxifloxacin"]
        )
        moxi_annos = [a for a in result if a.drug == "moxifloxacin"]
        assert len(moxi_annos) == 1

    def test_no_duplicate_drugs_when_all_sources_overlap(self, make_annotation):
        """When original, gene_associated, and GeneDatabase all share the same drug, only one annotation."""
        anno = make_annotation(drug="rifampicin", confidence="Assoc w R", source="WHO")
        processor = VariantProcessor()
        result = processor._expand_annotations_for_all_drugs(
            original_annotations=[anno],
            gene_associated_drugs=["rifampicin"],
            gene_id="Rv0667",  # GeneDatabase has ["rifampicin"]
        )
        result = list(result)
        assert len(result) == 1
        assert result[0].source == "WHO"
        assert result[0].confidence == "Assoc w R"


class TestCreateSyntheticAnnotations:
    """Tests for VariantProcessor._create_synthetic_annotations"""

    def test_empty_new_drugs_returns_empty(self, make_annotation):
        """Empty new_drugs should return empty list."""
        template = make_annotation(drug="rifampin")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs=set(),
        )
        assert result == []

    def test_cartesian_product(self, make_annotation):
        """N drugs x M templates should produce N*M synthetics."""
        t1 = make_annotation(drug="rifampin")
        t2 = make_annotation(drug="isoniazid", source="Other")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[t1, t2],
            new_drugs={"bedaquiline", "clofazimine"},
            confidence="No WHO annotation",
        )
        assert len(result) == 4  # 2 drugs x 2 templates

    def test_drug_field_always_overridden(self, make_annotation):
        """Drug field should always be overridden."""
        template = make_annotation(drug="rifampin")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs={"bedaquiline"},
        )
        assert len(result) == 1
        assert result[0].drug == "bedaquiline"

    def test_confidence_overridden_when_non_none(self, make_annotation):
        """Confidence should be overridden when non-None."""
        template = make_annotation(drug="rifampin", confidence="Assoc w R")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs={"bedaquiline"},
            confidence="No WHO annotation",
        )
        assert result[0].confidence == "No WHO annotation"

    def test_source_overridden_when_non_none(self, make_annotation):
        """Source should be overridden when non-None."""
        template = make_annotation(drug="rifampin", source="WHO")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs={"bedaquiline"},
            source="Mutation effect for given drug is not in TBDB",
        )
        assert result[0].source == "Mutation effect for given drug is not in TBDB"

    def test_comment_overridden_when_non_none(self, make_annotation):
        """Comment should be overridden when non-None."""
        template = make_annotation(drug="rifampin", comment="original comment")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs={"bedaquiline"},
            comment="new comment",
        )
        assert result[0].comment == "new comment"

    def test_none_values_preserve_template(self, make_annotation):
        """None override values should preserve template's original field value."""
        template = make_annotation(
            drug="rifampin", confidence="Assoc w R", source="WHO", comment="test comment",
        )
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs={"bedaquiline"},
            # confidence, source, comment all default to None
        )
        assert result[0].drug == "bedaquiline"  # always overridden
        assert result[0].confidence == "Assoc w R"  # preserved
        assert result[0].source == "WHO"  # preserved
        assert result[0].comment == "test comment"  # preserved

    def test_template_not_mutated(self, make_annotation):
        """Template objects should not be mutated."""
        template = make_annotation(drug="rifampin", confidence="Assoc w R")
        processor = VariantProcessor()
        processor._create_synthetic_annotations(
            base_annotations=[template],
            new_drugs={"bedaquiline"},
            confidence="No WHO annotation",
        )
        assert template.drug == "rifampin"
        assert template.confidence == "Assoc w R"

    def test_multiple_templates_per_drug(self, make_annotation):
        """Each drug should produce one synthetic per template."""
        t1 = make_annotation(drug="rifampin", source="WHO")
        t2 = make_annotation(drug="isoniazid", source="Expert")
        processor = VariantProcessor()
        result = processor._create_synthetic_annotations(
            base_annotations=[t1, t2],
            new_drugs={"bedaquiline"},
            confidence="No WHO annotation",
        )
        assert len(result) == 2
        assert all(r.drug == "bedaquiline" for r in result)
        sources = {r.source for r in result}
        assert "WHO" in sources
        assert "Expert" in sources


class TestDeduplicateVariants:
    """Tests for VariantProcessor.deduplicate_variants"""

    def test_no_duplicates_returns_all(self, make_variant):
        """No duplicates should return all variants unchanged."""
        v1 = make_variant(drug="rifampin")
        v2 = make_variant(drug="isoniazid")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v1, v2])
        assert len(result) == 2

    def test_empty_list_returns_empty(self):
        """Empty list should return empty list."""
        processor = VariantProcessor()
        result = processor.deduplicate_variants([])
        assert result == []

    @pytest.mark.parametrize("higher,lower", [
        ("Assoc w R", "Assoc w R - interim"),
        ("Assoc w R", "Uncertain significance"),
        ("Assoc w R - interim", "Not assoc w R - Interim"),
        ("Uncertain significance", "Not assoc w R"),
        ("Not assoc w R", "No WHO annotation"),
    ])
    def test_higher_confidence_replaces_lower(self, make_variant, higher, lower):
        """Higher confidence should replace lower for same variant."""
        v_high = make_variant(confidence=higher, source="WHO")
        v_low = make_variant(confidence=lower, source="WHO")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v_low, v_high])
        assert len(result) == 1
        assert result[0].confidence == higher

    @pytest.mark.parametrize("higher,lower", [
        ("Assoc w R", "Assoc w R - interim"),
        ("Assoc w R", "Uncertain significance"),
        ("Assoc w R - interim", "Not assoc w R - Interim"),
        ("Uncertain significance", "Not assoc w R"),
        ("Not assoc w R", "No WHO annotation"),
    ])
    def test_higher_confidence_kept_when_first(self, make_variant, higher, lower):
        """Higher confidence should be kept even when it appears first."""
        v_high = make_variant(confidence=higher, source="WHO")
        v_low = make_variant(confidence=lower, source="WHO")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v_high, v_low])
        assert len(result) == 1
        assert result[0].confidence == higher

    def test_who_source_preferred_at_equal_rank(self, make_variant):
        """WHO source should be preferred at equal confidence rank."""
        v_who = make_variant(confidence="Assoc w R", source="WHO")
        v_other = make_variant(confidence="Assoc w R", source="Expert")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v_other, v_who])
        assert len(result) == 1
        assert result[0].source == "WHO"

    def test_different_drugs_not_deduplicated(self, make_variant):
        """Different drugs should not be considered duplicates."""
        v1 = make_variant(drug="rifampin")
        v2 = make_variant(drug="isoniazid")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v1, v2])
        assert len(result) == 2

    def test_different_nucleotide_change_not_deduplicated(self, make_variant):
        """Different nucleotide_change should not be considered duplicates."""
        v1 = make_variant(nucleotide_change="c.1349C>T")
        v2 = make_variant(nucleotide_change="c.1350A>G")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v1, v2])
        assert len(result) == 2

    def test_winning_variant_retains_all_fields(self, make_variant):
        """Winning variant should retain all its original fields."""
        v_high = make_variant(confidence="Assoc w R", source="WHO", comment="high quality")
        v_low = make_variant(confidence="No WHO annotation", source="other")
        processor = VariantProcessor()
        result = processor.deduplicate_variants([v_low, v_high])
        assert len(result) == 1
        assert result[0].confidence == "Assoc w R"
        assert result[0].source == "WHO"
        assert result[0].comment == "high quality"


class TestGenerateUnreportedVariants:
    """Tests for VariantProcessor.generate_unreported_variants"""

    def test_empty_variants_generates_all(self):
        """Empty variants should generate all gene-drug combos from GeneDatabase."""
        processor = VariantProcessor()
        result = processor.generate_unreported_variants([], "test_sample")
        expected_count = sum(len(v["drugs"]) for v in GeneDatabase.GENE_DATABASE.values())
        assert len(result) == expected_count

    def test_covering_one_gene_reduces_count(self, make_variant):
        """Having a variant for one gene should reduce unreported count."""
        processor = VariantProcessor()
        v = make_variant(gene_id="Rv0667", gene_name="rpoB", drug="rifampin")
        all_unreported = processor.generate_unreported_variants([], "test_sample")
        partial_unreported = processor.generate_unreported_variants([v], "test_sample")
        rpoB_drug_count = len(GeneDatabase.get_drugs("Rv0667"))
        assert len(partial_unreported) == len(all_unreported) - rpoB_drug_count

    def test_all_genes_covered_returns_empty(self, make_variant):
        """If all GeneDatabase genes have variants, unreported should be empty."""
        processor = VariantProcessor()
        variants = []
        for gene_id in GeneDatabase.GENE_DATABASE:
            gene_name = GeneDatabase.get_gene_name(gene_id)
            drug = GeneDatabase.get_drugs(gene_id)[0]
            variants.append(make_variant(gene_id=gene_id, gene_name=gene_name, drug=drug))
        result = processor.generate_unreported_variants(variants, "test_sample")
        assert result == []

    def test_unreported_variants_have_thin_air_fields(self):
        """Unreported variants should have from_thin_air fields (pos='NA', etc.)."""
        processor = VariantProcessor()
        result = processor.generate_unreported_variants([], "test_sample")
        assert len(result) > 0
        for v in result:
            assert v.pos == "NA"
            assert v.depth == "NA"
            assert v.freq == "NA"
            assert v.type == "NA"

    def test_correct_sample_id_propagated(self):
        """Unreported variants should have the correct sample_id."""
        processor = VariantProcessor()
        result = processor.generate_unreported_variants([], "my_sample_123")
        assert all(v.sample_id == "my_sample_123" for v in result)

    def test_multi_drug_genes_generate_multiple_variants(self):
        """Genes with multiple drugs should generate multiple unreported variants."""
        processor = VariantProcessor()
        result = processor.generate_unreported_variants([], "test_sample")
        # gyrB (Rv0005) has 2 drugs: levofloxacin, moxifloxacin
        gyrB_unreported = [v for v in result if v.gene_id == "Rv0005"]
        assert len(gyrB_unreported) == 2

    def test_multiple_variants_same_gene_still_covered(self, make_variant):
        """Multiple variants for same gene_id should count as covered."""
        processor = VariantProcessor()
        v1 = make_variant(gene_id="Rv0667", gene_name="rpoB", drug="rifampin", nucleotide_change="c.2C>T")
        v2 = make_variant(gene_id="Rv0667", gene_name="rpoB", drug="rifampin", nucleotide_change="c.1350A>G")
        all_unreported = processor.generate_unreported_variants([], "test_sample")
        partial_unreported = processor.generate_unreported_variants([v1, v2], "test_sample")
        rpoB_drug_count = len(GeneDatabase.get_drugs("Rv0667"))
        assert len(partial_unreported) == len(all_unreported) - rpoB_drug_count


class TestProcessVariantRecords:
    """Tests for VariantProcessor.process_variant_records"""

    def test_empty_list_returns_empty(self):
        """Empty list should return empty list."""
        processor = VariantProcessor()
        result = processor.process_variant_records([])
        assert result == []

    def test_single_record_produces_variants(self, make_variant_record, make_annotation):
        """Single record with annotation should produce variants."""
        vr = make_variant_record(
            annotation=[make_annotation(drug="rifampin")],
            gene_associated_drugs=["rifampin"],
        )
        processor = VariantProcessor()
        result = processor.process_variant_records([vr])
        assert len(result) > 0
        drugs = {v.drug for v in result}
        assert "rifampin" in drugs

    def test_record_with_no_annotations_skipped(self, make_variant_record):
        """Record with no annotations should produce no variants."""
        vr = make_variant_record(annotation=[])
        processor = VariantProcessor()
        result = processor.process_variant_records([vr])
        assert result == []
