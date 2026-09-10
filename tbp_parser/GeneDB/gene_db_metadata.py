"""
The `--db_bed` file carries neither `tier` nor `promoter_region`, and it has no concept of locus tag
aliases, so `build_gene_db` cannot derive any of them from its input currently. Genes absent from this map
default to a tier of "NA", an empty promoter region, and no aliases.

Sources
----------------
Both fields originate from the WHO catalogue of mutations in M. tuberculosis, 2nd edition (2023):

- `tier`            WHO v2 catalogue, Table 21 "Candidate resistance genes" (p. 89)
- `promoter_region` WHO v2 catalogue, Table 22 "Upstream/promoter regions of candidate resistance genes" (p. 89, 90)

The same data is published as CSVs in WHO's companion repository (https://github.com/GTB-tbsequencing/mutation-catalogue-2023).
See `Input data files for Solo algorithms/additional-data/gene_promoters.csv` (matches Table 22).
Tiers are not published as a flat file there, so Table 21 and the `tier` column of
`Final Result Files/WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt` are the only public sources.

Promoter Region Coordinates
----------------
Table 22 lists each region as `1-N` (e.g. aftB `1-129`), counted back from the gene start (CDS 5' end)
found in `gene_locations.csv`. Those are stored here as negative ranges (`[-129, -1]`), matching HGVS
`c.-N` numbering back from the start codon.

Table 22's "primary transcriptional start site" column is what set each region's length. It says,
"for each resistance gene, relevant promoter and/or upstream regions were defined according to the primary
transcriptional start site" (p. 88), with the region truncated if it extended into an adjacent coding sequence.

The offset of bp derived from Table 22 is exactly 51 bp: `region_end - (gene_start - TSS)` Genes with no listed
TSS fall back to a bare 51 bp (`[-51, -1]` entries).

Tier Definition
----------------
A tier is WHO's predefined ranking of where resistance mutations are expected to occur for each drug.
Tier 1 covers the gene sequences and promoters "considered most likely to contain resistance mutations".
Tier 2 covers the remaining candidate genes, with "a lower, but still reasonable pre-test probability".
Only these two tiers exist.

Caveats
----------------
1. Tiers are technically assigned per (drug, gene), not per gene: 9 genes hold both tiers depending on the drug
   (Rv1258c, Rv1979c, Rv2983, fbiA, fbiB, fbiC, fgd1, mshA, rrl), and a gene may be a candidate for
   one drug and untiered for another.
2. `Rv0678` is absent from Table 22: WHO folds its upstream region into mmpS5's `1-85` window
   (Table 22 footnote d). The `[-84, -1]` entry here is tbp-parser splitting that shared intergenic
   region back out. Likewise aftA, fabG1 and furA have no entries of their own because Table 22
   subsumes them into the embC, inhA and katG promoter regions respectively (footnotes a, b, c).

"""

GENE_DB_METADATA: dict[str, dict] = {
    "Rv0001": {"tier": "NA", "promoter_region": [-314, -1]},  # dnaA
    "Rv0005": {"tier": "Tier 1", "promoter_region": [-108, -1]},  # gyrB
    "Rv0006": {"tier": "Tier 1", "promoter_region": [-35, -1]},  # gyrA
    "Rv0010c": {"tier": "NA", "promoter_region": [-156, -1]},  # Rv0010c
    "Rv0407": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # fgd1
    "Rv0486": {"tier": "Tier 2", "promoter_region": [-669, -1]},  # mshA
    "Rv0529": {"tier": "Tier 2", "promoter_region": [-191, -1]},  # ccsA
    "Rv0565c": {"tier": "NA", "promoter_region": [-78, -1]},  # Rv0565c
    "Rv0635": {"tier": "NA", "promoter_region": [-51, -1]},  # hadA
    "Rv0639": {"tier": "NA", "promoter_region": [-201, -1]},  # nusG
    "Rv0643c": {"tier": "NA", "promoter_region": []},  # mmaA3
    "Rv0667": {"tier": "Tier 1", "promoter_region": [-263, -1]},  # rpoB
    "Rv0668": {"tier": "Tier 2", "promoter_region": [-45, -1]},  # rpoC
    "Rv0676c": {"tier": "Tier 1", "promoter_region": []},  # mmpL5
    "Rv0677c": {"tier": "Tier 1", "promoter_region": [-85, -1]},  # mmpS5
    "Rv0678": {"tier": "Tier 1", "promoter_region": [-84, -1]},  # mmpR5
    "Rv0682": {"tier": "Tier 1", "promoter_region": [-234, -1]},  # rpsL
    "Rv0701": {"tier": "Tier 1", "promoter_region": [[-51, -1], [-503, -323]]},  # rplC
    "Rv1129c": {"tier": "NA", "promoter_region": [-51, -1]},  # Rv1129c
    "Rv1173": {"tier": "Tier 1", "promoter_region": [-127, -1]},  # fbiC
    "Rv1221": {"tier": "NA", "promoter_region": [-51, -1]},  # sigE
    "Rv1258c": {"tier": "Tier 1", "promoter_region": [-58, -1]},  # Rv1258c
    "Rv1267c": {"tier": "Tier 2", "promoter_region": [-103, -1]},  # embR
    "Rv1305": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # atpE
    "EBG00000313325": {"tier": "Tier 1", "promoter_region": [-151, -1], "aliases": ["Rvnr01", "MTB000019"]},  # rrs
    "EBG00000313339": {"tier": "Tier 1", "promoter_region": [-51, -1], "aliases": ["Rvnr02", "MTB000020"]},  # rrl
    "Rv1484": {"tier": "Tier 1", "promoter_region": [-813, -1]},  # inhA
    "Rv1630": {"tier": "NA", "promoter_region": [-100, -1]},  # rpsA
    "Rv1644": {"tier": "NA", "promoter_region": [-51, -1]},  # tsnR
    "Rv1694": {"tier": "Tier 1", "promoter_region": [[-51, -1], [-236, -185]]},  # tlyA
    "Rv1819c": {"tier": "NA", "promoter_region": [-81, -1]},  # bacA
    "Rv1854c": {"tier": "Tier 2", "promoter_region": [-96, -1]},  # ndh
    "Rv1908c": {"tier": "Tier 1", "promoter_region": [-532, -1]},  # katG
    "Rv1918c": {"tier": "Tier 2", "promoter_region": [-122, -1]},  # PPE35
    "Rv1979c": {"tier": "Tier 2", "promoter_region": [-470, -1]},  # Rv1979c
    "Rv2043c": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # pncA
    "Rv2245": {"tier": "NA", "promoter_region": []},  # kasA
    "Rv2416c": {"tier": "Tier 1", "promoter_region": [-84, -1]},  # eis
    "Rv2428": {"tier": "Tier 1", "promoter_region": [-93, -1]},  # ahpC
    "Rv2447c": {"tier": "NA", "promoter_region": []},  # folC
    "Rv2477c": {"tier": "NA", "promoter_region": [-88, -1]},  # Rv2477c
    "Rv2535c": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # pepQ
    "Rv2671": {"tier": "NA", "promoter_region": []},  # ribD
    "Rv2680": {"tier": "NA", "promoter_region": [-153, -1]},  # Rv2680
    "Rv2681": {"tier": "NA", "promoter_region": [-2, -1]},  # Rv2681
    "Rv2752c": {"tier": "Tier 2", "promoter_region": [[-51, -1], [-984, -934]]},  # Rv2752c
    "Rv2754c": {"tier": "NA", "promoter_region": []},  # thyX
    "Rv2764c": {"tier": "NA", "promoter_region": []},  # thyA
    "Rv2780": {"tier": "NA", "promoter_region": []},  # ald
    "Rv2983": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # fbiD
    "Rv3083": {"tier": "Tier 2", "promoter_region": [-51, -1]},  # Rv3083
    "Rv3197A": {"tier": "Tier 1", "promoter_region": [-404, -1]},  # whiB7
    "Rv3236c": {"tier": "Tier 2", "promoter_region": [[-51, -1], [-538, -488]]},  # Rv3236c
    "Rv3244c": {"tier": "NA", "promoter_region": []},  # lpqB
    "Rv3245c": {"tier": "NA", "promoter_region": [-50, -1]},  # mtrB
    "Rv3246c": {"tier": "NA", "promoter_region": [-376, -1]},  # mtrA
    "Rv3261": {"tier": "Tier 1", "promoter_region": [-138, -1]},  # fbiA
    "Rv3262": {"tier": "Tier 1", "promoter_region": []},  # fbiB
    "Rv3423c": {"tier": "NA", "promoter_region": []},  # alr
    "Rv3457c": {"tier": "Tier 2", "promoter_region": [-536, -1]},  # rpoA
    "Rv3547": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # ddn
    "Rv3596c": {"tier": "Tier 1", "promoter_region": [-106, -1]},  # clpC1
    "Rv3601c": {"tier": "Tier 1", "promoter_region": [[-51, -1], [-1949, -1838]]},  # panD
    "Rv3696c": {"tier": "NA", "promoter_region": [-52, -1]},  # glpK
    "Rv3793": {"tier": "Tier 1", "promoter_region": [-1982, -1]},  # embC
    "Rv3794": {"tier": "Tier 1", "promoter_region": [-86, -1]},  # embA
    "Rv3795": {"tier": "Tier 1", "promoter_region": []},  # embB
    "Rv3805c": {"tier": "Tier 2", "promoter_region": [-129, -1]},  # aftB
    "Rv3806c": {"tier": "Tier 2", "promoter_region": [-51, -1]},  # ubiA
    "Rv3854c": {"tier": "Tier 1", "promoter_region": [-51, -1]},  # ethA
    "Rv3855": {"tier": "Tier 2", "promoter_region": [-26, -1]},  # ethR
    "Rv3862c": {"tier": "Tier 2", "promoter_region": [-126, -1]},  # whiB6
    "Rv3919c": {"tier": "Tier 1", "promoter_region": [-79, -1]},  # gid
}
