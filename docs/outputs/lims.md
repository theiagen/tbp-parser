---
title: LIMS Report
---

The LIMS report is intended for direct import into a STAR LIMS system. The columns are in the specific LIMS code format for CDPH, and may not apply to your LIMS system. Please contact us if you need different column headers and we can work with you towards a solution.

### Explanation of column headers

| Column name | Explanation |
| --- | --- |
| MDL sample accession numbers | The name of the sample |
| M_DST_A01_ID | The lineage of the sample in human-readable language |
| M_DST_B01_INH | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (isoniazid) |
| M_DST_B02_katG | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for isoniazid |
| M_DST_B03_fabG1 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for isoniazid|
| M_DST_B04_inhA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for isoniazid |
| M_DST_C01_ETO | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (ethionamide) |
| M_DST_C02_ethA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for ethionamide |
| M_DST_C03_fabG1 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for ethionamide |
| M_DST_C04_inhA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for ethionamide |
| M_DST_D01_RIF | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (rifampin) |
| M_DST_D02_rpoB | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for rifampin |
| M_DST_E01_PZA | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (pyrazinamide) |
| M_DST_E02_pncA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for pyrazinamide |
| M_DST_F01_EMB | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (ethambutol) |
| M_DST_F02_embA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for ethambutol |
| M_DST_F03_embB | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for ethambutol |
| M_DST_G01_AMK | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (amikacin) |
| M_DST_G02_rrs | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for amikacin |
| M_DST_G03_eis | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for amikacin |
| M_DST_H01_KAN | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (kanamycin) |
| M_DST_H02_rrs | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for kanamycin |
| M_DST_H03_eis | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for kanamycin |
| M_DST_I01_CAP | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (capreomycin) |
| M_DST_I02_rrs | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for capreomycin |
| M_DST_I03_tlyA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for capreomycin |
| M_DST_J01_MFX | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (moxifloxacin) |
| M_DST_J02_gyrA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for moxifloxacin |
| M_DST_J03_gyrB | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for moxifloxacin |
| M_DST_K01_LFX | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (levofloxacin) |
| M_DST_K02_gyrA | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for levofloxacin |
| M_DST_K03_gyrB | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for levofloxacin |
| M_DST_L01_BDQ | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (bedaquiline) |
| M_DST_L02_Rv0678 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for bedaquiline |
| M_DST_L03_atpE | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for bedaquiline |
| M_DST_L04_pepQ | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for bedaquiline |
| M_DST_L05_mmpL5 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for bedaquiline |
| M_DST_L06_mmpS5 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for bedaquiline |
| M_DST_M01_CFZ | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (clofazimine) |
| M_DST_M02_Rv0678 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for clofazimine |
| M_DST_M03_pepQ | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for clofazimine |
| M_DST_M04_mmpL5 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for clofazimine |
| M_DST_M05_mmpS5 | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for clofazimine |
| M_DST_N01_LZD | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (linezolid) |
| M_DST_N02_rrl | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for linezolid |
| M_DST_N03_rplC | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for linezolid |
| Analysis date | The date `tbp-parser` was run in YYYY-MM-DD HH:SS format |
| Operator | The name of the person who ran `tbp-parser`; can be provided with the `--operator` input parameter. If left blank, “Operator not provided” is the default value. |
| M_DST_O01_lineage | The lineage of the sample (the `main_lin` of the sample as reported by TBProfiler) |
| M_DST_P01_CS | The highest `mdl_interpretation` resistance identified for mutations associated with this drug (cycloserine); only included when `--add_cs_lims` is set to true |
| M_DST_P02_ald | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for cycloserine; only included when `--add_cs_lims` is set to true |
| M_DST_PO3_alr | Any non-S mutations found in this gene with good quality responsible for the predicted resistance for cycloserine; only included when `--add_cs_lims` is set to true |

The LIMS report offers a condensed version of the laboratorian report with more details than the Looker report. By containing only the most important information about a drug and its related mutations, the LIMS report provides an invaluable summary.
