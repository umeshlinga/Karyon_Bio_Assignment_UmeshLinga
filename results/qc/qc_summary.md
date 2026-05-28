# GSE136103 — QC Summary

- Samples discovered: **20**
- Cells loaded (post-subsample): **28,173**
- Cells retained after QC: **26,797** (95.1%)
- Predicted doublets (Scrublet, pre-filter): **201**

## QC thresholds (liver-fibrosis-aware)
- min_genes >= 200; min_counts >= 500
- pct_counts_mt <= 20.0%  *(relaxed: activated HSCs and SAMs in fibrotic livers legitimately have elevated mito)*
- pct_counts_hb <= 50.0%; pct_counts_ribo <= 60.0%

## Per-sample retention

| sample_id                    |   n_cells_pre |   med_genes |   med_counts |   med_mt |   n_cells_post |   retention_pct |
|:-----------------------------|--------------:|------------:|-------------:|---------:|---------------:|----------------:|
| GSM4041150_healthy1_cd45+    |          1500 |      1038.5 |       2680   |      3.2 |           1480 |            98.7 |
| GSM4041151_healthy1_cd45-A   |          1039 |      1204   |       3485   |      5.9 |            816 |            78.5 |
| GSM4041152_healthy1_cd45-B   |           478 |      1369.5 |       3381.5 |      6.9 |            371 |            77.6 |
| GSM4041153_healthy2_cd45+    |          1500 |       588   |       1128.5 |      2.3 |           1497 |            99.8 |
| GSM4041154_healthy2_cd45-    |          1500 |       677.5 |       1359   |      2.7 |           1486 |            99.1 |
| GSM4041155_healthy3_cd45+    |          1500 |      1352   |       3953   |      3.8 |           1462 |            97.5 |
| GSM4041156_healthy3_cd45-A   |          1500 |      1756   |       4694   |      4.2 |           1397 |            93.1 |
| GSM4041157_healthy3_cd45-B   |          1156 |      1817   |       5542   |      4.4 |           1074 |            92.9 |
| GSM4041158_healthy4_cd45+    |          1500 |      1234.5 |       3814   |      3.6 |           1481 |            98.7 |
| GSM4041159_healthy4_cd45-    |          1500 |      1927.5 |       6269.5 |      3.7 |           1470 |            98   |
| GSM4041160_healthy5_cd45+    |          1500 |       905.5 |       2262.5 |      2.3 |           1490 |            99.3 |
| GSM4041161_cirrhotic1_cd45+  |          1500 |      1154   |       3475.5 |      4.1 |           1402 |            93.5 |
| GSM4041162_cirrhotic1_cd45-A |          1500 |      1607   |       4838.5 |      5.2 |           1367 |            91.1 |
| GSM4041163_cirrhotic1_cd45-B |          1500 |      1538   |       4854   |      6.3 |           1302 |            86.8 |
| GSM4041164_cirrhotic2_cd45+  |          1500 |      1333   |       3319.5 |      3.4 |           1444 |            96.3 |
| GSM4041165_cirrhotic2_cd45-  |          1500 |      1744   |       4852.5 |      3.8 |           1429 |            95.3 |
| GSM4041166_cirrhotic3_cd45+  |          1500 |      1297   |       4057   |      4.3 |           1401 |            93.4 |
| GSM4041167_cirrhotic3_cd45-  |          1500 |      1463.5 |       4457   |      4.2 |           1443 |            96.2 |
| GSM4041168_cirrhotic4_cd45+  |          1500 |      1102.5 |       2937   |      2.5 |           1495 |            99.7 |
| GSM4041169_cirrhotic5_cd45+  |          1500 |      1524   |       5330   |      3.2 |           1490 |            99.3 |
