# M8 HERB Ingredient-Target Acquisition Run

## Scope

This run downloads HERB ingredient detail payloads through the public JSON detail API and extracts compound-target relationships. It separates literature-mined PubMed sentence-supported targets from database-integrated ingredient targets.

## API

- Endpoint: `http://47.92.70.12/chedi/api/?`
- Method: POST JSON body with `func_name=detail_api`, `label=Ingredient`, `v=<Ingredient_id>`, and `key_id=<Ingredient_id>`.

## Outputs

- `data/raw/herb2/ingredient_details/*.json`
- `results/m8_herb_ingredient_targets/herb_ingredient_target_edges.tsv`
- `results/m8_herb_ingredient_targets/herb_ingredient_target_summary.tsv`
- `results/m8_herb_ingredient_targets/herb_detail_download_audit.tsv`

## Acquisition Summary

- Requested ingredients: 27
- Successful detail payloads: 27
- Failed detail payloads: 0
- Literature/PubMed target edge rows: 3767
- Database-integrated target edge rows: 3829

## Top Candidates By Disease-Module Target Overlap

| Ingredient | Literature targets | Database targets | UC/IBD up150 target hits | UC/IBD down150 target hits | Up-hit genes |
|---|---:|---:|---:|---:|---|
| Trans-resveratrol | 214 | 1513 | 52 | 15 | ADGRE2;ARNTL2;BCL2A1;C3;C4BPB;CCL2;CHI3L1;CLDN1;COL1A1;CXCL1;CXCL10;CXCL6;CXCL8;DEFA5;DUOX2;FOS;FPR1;ICAM1;IDO1;IGFBP5;IL1A;IL1B;IL1RN;IL7R;ITGA5;MMP1;MMP3;MMP7;MMP9;NOS2;PDZK1IP1;PLAU;PTGS2;REG3A;RGS5;S100P;SELE;SERPINA1;SERPINA3;SERPINB3;SERPINB5;SERPINB7;SERPINB9;SLC2A3;SLC7A11;SLCO1B3;SOCS3;SPP1;TDO2;TGM2;TIMP1;TNFAIP6 |
| Curcumin | 307 | 152 | 23 | 4 | CCL2;CD274;COL1A1;CTLA4;CXCL1;CXCL8;FOS;ICAM1;IL1A;IL1B;IL33;MMP1;MMP3;MMP9;NOS2;PLAU;PTGS2;S100A8;SELE;SELP;SERPINA1;SPP1;TIMP1 |
| Quercetin | 189 | 161 | 23 | 5 | CCL2;COL1A1;CXCL10;CXCL11;CXCL8;DUOX2;FOS;ICAM1;IDO1;IL1A;IL1B;MMP1;MMP12;MMP3;MMP9;NOS2;PLAU;PTGS2;SELE;SLC7A11;SPP1;STC1;TREM1 |
| Wogonin | 73 | 308 | 16 | 4 | BCL2A1;C2;CCL2;CFB;CFI;CXCL8;GZMK;IL1B;KLK10;MMP1;MMP3;MMP9;NOS2;PLAU;PTGS2;SPP1 |
| Capsaicin | 60 | 334 | 15 | 4 | C2;CFB;CFI;CXCL8;FOS;GZMK;ICAM1;IL1B;KLK10;MMP1;MMP9;NOS2;PLAU;PTGS2;TIMP1 |
| Melatonin | 161 | 106 | 14 | 3 | CD274;COL1A1;CXCL8;FOS;ICAM1;IL1A;IL1B;LCN2;MMP3;MMP9;NNMT;NOS2;PTGS2;SLCO1B3 |
| Archin | 96 | 263 | 14 | 5 | BCL2A1;CCL2;CXCL1;CXCL5;CXCL8;IL1A;IL1B;MMP1;MMP3;MMP9;NOS2;PLAU;PTGS2;SLCO1B3 |
| Berberine | 187 | 82 | 13 | 7 | CCL2;CD274;CXCL8;FOS;IL1A;IL1B;MMP3;MMP9;NOS2;OSMR;PTGS2;TIMP1;VWF |
| Berberime | 77 | 334 | 13 | 3 | FAP;FOS;ICAM1;IDO1;IL1B;MMP1;MMP3;MMP9;NOS2;PLAU;PTGS2;SLCO1B3;TIMP1 |
| Polydatin | 54 | 18 | 10 | 2 | CCL2;CXCL10;CXCL8;ICAM1;IL1B;NOS2;PTGS2;SELE;SELP;SPP1 |
| Cis-resveratrol | 152 | 37 | 9 | 4 | CXCL10;FOS;ICAM1;IL1B;MMP10;MMP9;NOS2;PTGS2;SELP |
| Ellagic acid | 77 | 18 | 8 | 1 | CCL2;CTSK;CXCL8;FOS;IL1B;MMP9;NOS2;PTGS2 |

## Interpretation Boundary

- `drug_paper_target` is the preferred HERB target layer because it carries PubMed IDs, relationship labels, grades, and supporting sentences.
- `ingredient_target` is retained as a secondary annotation layer because it is database-integrated and may include predicted or inherited targets.
- Target overlap with the bulk UC/IBD module is hypothesis-generating. It should be combined with scRNA localization, GWAS/MR or independent disease-gene evidence before mechanism claims.
