# Full CellChatDB.mouse interpretation (GSE274279, 2,019 interactions)

## Communication strength by age

| age | nonzero edges | significant edges | significant interactions | significant pathways | total significant probability |
| --- | --- | --- | --- | --- | --- |
| 3M | 739 | 455 | 146 | 33 | 11.49 |
| 12M | 561 | 360 | 115 | 34 | 9.392 |
| 24M | 463 | 245 | 64 | 27 | 7.445 |

## Persistence of significant interactions

- significant in all three ages: 46
- significant in 3M only: 0
- significant in 24M only: 2
- significant in 3M and 12M: 0

Core interactions significant at every age (46):
  AGRN_DAG1, BMP5_ACVR1_ACVR2A, BMP5_ACVR1_BMPR2, BMP5_BMPR1A_ACVR2A, BMP5_BMPR1A_BMPR2, BMP5_BMPR1B_ACVR2A, BMP5_BMPR1B_BMPR2, BMP6_ACVR1_ACVR2A, BMP6_BMPR1A_ACVR2A, BMP6_BMPR1B_ACVR2A, BMP6_BMPR1B_BMPR2, CADM1_CADM1, COL4A3_ITGA9_ITGB1, COL4A5_ITGA9_ITGB1, EFNA5_EPHA3, EFNA5_EPHA4, EFNA5_EPHA5, EFNA5_EPHA7, EFNA5_EPHB2, FN1_ITGA8_ITGB1, GRN_SORT1, LAMA2_ITGA9_ITGB1, LAMA4_DAG1, LAMA4_ITGA9_ITGB1, LAMB1_DAG1, LAMB1_ITGA9_ITGB1, LAMC1_ITGA9_ITGB1, LRRC4C_NTNG1, MPZL1_MPZL1, NCAM1_FGFR1, NCAM1_NCAM1, NEGR1_NEGR1, NRG3_ERBB4, NRXN1_NLGN1, NRXN3_NLGN1, PTN_ALK, PTN_PTPRZ1, PTN_SDC2, PTPRC_MRC1, PTPRM_PTPRM, SEMA3A_NRP1_PLXNA2, SEMA3A_NRP1_PLXNA4, SPP1_ITGA8_ITGB1, SPP1_ITGA9_ITGB1, TGFB2_ACVR1_TGFBR1, TGFB2_TGFBR1_TGFBR2

## Pathways with the largest total probability per age

### 3M
- PTN: 3.465 (17 edges)
- EPHA: 2.065 (20 edges)
- NRG: 1.692 (7 edges)
- PTPRM: 0.5261 (6 edges)
- LAMININ: 0.4867 (17 edges)
- NRXN: 0.4578 (10 edges)
- BMP: 0.4141 (15 edges)
- FGF: 0.3357 (7 edges)

### 12M
- PTN: 4.265 (29 edges)
- EPHA: 1.485 (17 edges)
- NRG: 0.9884 (7 edges)
- APP: 0.9688 (7 edges)
- BMP: 0.3239 (17 edges)
- NRXN: 0.2173 (10 edges)
- LAMININ: 0.1824 (21 edges)
- FN1: 0.1081 (8 edges)

### 24M
- PTN: 2.446 (16 edges)
- EPHA: 1.489 (16 edges)
- NRG: 1.384 (7 edges)
- NRXN: 0.4566 (12 edges)
- PTPRM: 0.4199 (12 edges)
- APP: 0.2196 (7 edges)
- BMP: 0.1825 (22 edges)
- SEMA3: 0.1421 (6 edges)

## Communication strength by cell type (significant edges)

### 3M
- Spiral ganglion neurons: outgoing 3.636 | incoming 2.359
- Glia/Schwann: outgoing 2.006 | incoming 2.38
- Supporting cells: outgoing 1.527 | incoming 1.978
- Fibrocytes: outgoing 1.509 | incoming 0.8305
- Inner hair cells: outgoing 1.243 | incoming 3.024
- Outer hair cells: outgoing 1.023 | incoming 0.17
- Macrophages/Microglia: outgoing 0.5474 | incoming 0.7495

### 12M
- Spiral ganglion neurons: outgoing 2.036 | incoming 0.3633
- Glia/Schwann: outgoing 1.885 | incoming 2.702
- Fibrocytes: outgoing 1.496 | incoming 1.022
- Supporting cells: outgoing 1.333 | incoming 0.9975
- Outer hair cells: outgoing 0.9166 | incoming 0.2955
- Inner hair cells: outgoing 0.8836 | incoming 2.264
- Macrophages/Microglia: outgoing 0.8424 | incoming 1.748

### 24M
- Spiral ganglion neurons: outgoing 2.169 | incoming 0.6929
- Glia/Schwann: outgoing 1.526 | incoming 2.039
- Supporting cells: outgoing 0.9653 | incoming 0.8493
- Fibrocytes: outgoing 0.8439 | incoming 0.6971
- Inner hair cells: outgoing 0.8345 | incoming 2.012
- Outer hair cells: outgoing 0.5855 | incoming 0.4762
- Macrophages/Microglia: outgoing 0.5208 | incoming 0.6786

## Focus interactions (manuscript candidates)

- AGRN_DAG1: 3M sig (prob 0.00151, p 0.01); 3M sig (prob 0.00144, p 0); 12M sig (prob 0.00277, p 0); 12M sig (prob 0.00275, p 0); 12M sig (prob 0.00161, p 0); 24M sig (prob 0.00201, p 0); 24M sig (prob 0.000687, p 0)
- APP_CD74: 12M sig (prob 0.152, p 0); 12M sig (prob 0.151, p 0); 12M sig (prob 0.148, p 0); 12M sig (prob 0.144, p 0); 12M sig (prob 0.138, p 0); 12M sig (prob 0.122, p 0); 12M sig (prob 0.114, p 0); 24M sig (prob 0.0384, p 0); 24M sig (prob 0.0381, p 0); 24M sig (prob 0.0354, p 0); 24M sig (prob 0.0316, p 0); 24M sig (prob 0.0263, p 0); 24M sig (prob 0.0261, p 0); 24M sig (prob 0.0238, p 0)
- BMP5_ACVR1_ACVR2A: 3M sig (prob 0.00556, p 0); 3M sig (prob 0.00466, p 0); 3M sig (prob 0.00448, p 0.03); 3M sig (prob 0.00373, p 0); 3M sig (prob 0.00312, p 0); 3M ns (prob 0.003, p 0.12); 12M sig (prob 0.00491, p 0); 12M sig (prob 0.00445, p 0); 12M ns (prob 0.00267, p 0.11); 12M sig (prob 0.00242, p 0); 24M sig (prob 0.00524, p 0); 24M sig (prob 0.00424, p 0); 24M sig (prob 0.00369, p 0.04); 24M sig (prob 0.00362, p 0.01); 24M ns (prob 0.00257, p 0.2); 24M ns (prob 0.00208, p 0.24); 24M ns (prob 0.00181, p 0.19); 24M ns (prob 0.00177, p 0.24)
- CADM1_CADM1: 3M sig (prob 0.173, p 0); 3M sig (prob 0.0372, p 0); 3M sig (prob 0.0372, p 0); 3M sig (prob 0.0286, p 0); 3M sig (prob 0.0286, p 0); 3M sig (prob 0.00709, p 0.03); 3M ns (prob 0.00541, p 0.05); 3M ns (prob 0.00541, p 0.05); 3M ns (prob 0.00412, p 0.89); 12M sig (prob 0.00988, p 0); 12M sig (prob 0.00621, p 0); 12M sig (prob 0.00621, p 0); 12M sig (prob 0.00516, p 0); 12M sig (prob 0.00516, p 0); 12M sig (prob 0.00455, p 0); 12M sig (prob 0.00455, p 0); 12M sig (prob 0.00389, p 0); 12M sig (prob 0.00323, p 0); 12M sig (prob 0.00323, p 0); 12M sig (prob 0.00285, p 0); 12M sig (prob 0.00285, p 0); 12M ns (prob 0.00268, p 0.12); 12M sig (prob 0.00237, p 0); 12M sig (prob 0.00237, p 0); 12M sig (prob 0.00209, p 0); 24M sig (prob 0.0507, p 0); 24M sig (prob 0.0121, p 0); 24M sig (prob 0.0121, p 0); 24M ns (prob 0.00281, p 0.07)
- CDH2_CDH2: 3M sig (prob 0.0556, p 0); 3M sig (prob 0.0108, p 0); 3M sig (prob 0.0108, p 0); 3M sig (prob 0.00201, p 0); 24M sig (prob 0.00532, p 0)
- COL4A3_ITGA3_ITGB1: 3M sig (prob 0.00192, p 0)
- EFNA5_EPHA3: 3M sig (prob 0.165, p 0); 3M sig (prob 0.086, p 0); 3M sig (prob 0.0626, p 0); 3M sig (prob 0.0471, p 0); 3M sig (prob 0.0437, p 0.04); 3M sig (prob 0.0403, p 0.01); 3M ns (prob 0.0375, p 0.07); 3M ns (prob 0.0354, p 0.14); 3M ns (prob 0.0333, p 0.36); 3M ns (prob 0.029, p 0.48); 3M ns (prob 0.0263, p 1); 3M ns (prob 0.0256, p 0.58); 3M ns (prob 0.023, p 0.86); 3M ns (prob 0.0213, p 0.9); 3M ns (prob 0.0196, p 0.93); 3M ns (prob 0.0182, p 0.94); 3M ns (prob 0.0172, p 0.95); 3M ns (prob 0.0164, p 0.86); 3M ns (prob 0.0161, p 0.95); 3M ns (prob 0.0152, p 0.87); 3M ns (prob 0.014, p 0.88); 3M ns (prob 0.013, p 0.89); 3M ns (prob 0.0123, p 0.89); 3M ns (prob 0.0115, p 0.89); 3M ns (prob 0.00741, p 0.9); 3M ns (prob 0.00685, p 0.91); 3M ns (prob 0.0067, p 1); 3M ns (prob 0.00654, p 0.82); 3M ns (prob 0.0063, p 0.91); 3M ns (prob 0.0062, p 1); 3M ns (prob 0.00605, p 0.82); 3M ns (prob 0.00585, p 0.91); 3M ns (prob 0.0057, p 1); 3M ns (prob 0.00556, p 0.82); 3M ns (prob 0.00551, p 0.91); 3M ns (prob 0.00529, p 1); 3M ns (prob 0.00518, p 0.91); 3M ns (prob 0.00516, p 0.82); 3M ns (prob 0.00499, p 1); 3M ns (prob 0.00486, p 0.82); 3M ns (prob 0.00468, p 1); 3M ns (prob 0.00457, p 0.82); 12M sig (prob 0.129, p 0); 12M sig (prob 0.0559, p 0); 12M sig (prob 0.0399, p 0); 12M sig (prob 0.0391, p 0); 12M sig (prob 0.0339, p 0); 12M sig (prob 0.0336, p 0); 12M sig (prob 0.0322, p 0); 12M sig (prob 0.027, p 0); 12M sig (prob 0.0248, p 0.04); 12M ns (prob 0.0163, p 0.05); 12M sig (prob 0.0159, p 0.04); 12M ns (prob 0.0138, p 0.1); 12M ns (prob 0.0137, p 0.11); 12M ns (prob 0.0131, p 0.1); 12M ns (prob 0.0109, p 0.27); 12M ns (prob 0.00705, p 0.37); 12M ns (prob 0.00691, p 0.36); 12M ns (prob 0.00596, p 0.38); 12M ns (prob 0.00591, p 0.4); 12M ns (prob 0.00566, p 0.38); 12M ns (prob 0.00473, p 0.42); 24M sig (prob 0.137, p 0); 24M sig (prob 0.0445, p 0); 24M sig (prob 0.0415, p 0); 24M sig (prob 0.0407, p 0.02); 24M ns (prob 0.0375, p 0.06); 24M sig (prob 0.0349, p 0); 24M sig (prob 0.0346, p 0); 24M sig (prob 0.0267, p 0); 24M sig (prob 0.0247, p 0.02); 24M ns (prob 0.0123, p 0.27); 24M ns (prob 0.0114, p 0.25); 24M ns (prob 0.0113, p 0.29); 24M ns (prob 0.0105, p 0.3); 24M ns (prob 0.00953, p 0.31); 24M ns (prob 0.00944, p 0.32); 24M ns (prob 0.00878, p 0.33); 24M ns (prob 0.0087, p 0.33); 24M ns (prob 0.00726, p 0.3); 24M ns (prob 0.00671, p 0.33); 24M ns (prob 0.00669, p 0.36); 24M ns (prob 0.00618, p 0.34)
- EFNA5_EPHA4: 3M sig (prob 0.0482, p 0); 3M sig (prob 0.0255, p 0); 3M sig (prob 0.0243, p 0); 3M sig (prob 0.0235, p 0); 3M sig (prob 0.0168, p 0); 3M sig (prob 0.0151, p 0); 3M sig (prob 0.0123, p 0); 3M sig (prob 0.0117, p 0); 3M sig (prob 0.00878, p 0); 3M sig (prob 0.00833, p 0); 3M sig (prob 0.00758, p 0.02); 3M sig (prob 0.00724, p 0); 3M sig (prob 0.00686, p 0); 3M ns (prob 0.00669, p 0.12); 3M ns (prob 0.00515, p 0.39); 3M ns (prob 0.00394, p 0.62); 3M ns (prob 0.00374, p 0.69); 3M ns (prob 0.00356, p 0.93); 3M ns (prob 0.00348, p 0.61); 3M ns (prob 0.00338, p 1); 3M ns (prob 0.0033, p 0.69); 3M ns (prob 0.00231, p 0.86); 3M ns (prob 0.00209, p 1); 3M ns (prob 0.00204, p 0.81); 12M sig (prob 0.0221, p 0); 12M sig (prob 0.0184, p 0); 12M sig (prob 0.015, p 0); 12M sig (prob 0.00894, p 0); 12M sig (prob 0.00741, p 0); 12M sig (prob 0.00605, p 0); 12M ns (prob 0.00386, p 0.06); 12M ns (prob 0.0032, p 0.11); 12M ns (prob 0.00261, p 0.2); 24M sig (prob 0.022, p 0); 24M sig (prob 0.0175, p 0); 24M sig (prob 0.0158, p 0); 24M sig (prob 0.00596, p 0.02); 24M ns (prob 0.00549, p 0.05); 24M sig (prob 0.00471, p 0.03); 24M ns (prob 0.00434, p 0.09); 24M ns (prob 0.00424, p 0.05); 24M ns (prob 0.00391, p 0.16)
- EFNB2_EPHA4: 3M sig (prob 0.00634, p 0); 3M sig (prob 0.00329, p 0); 3M sig (prob 0.00312, p 0); 3M sig (prob 0.00193, p 0); 12M sig (prob 0.00246, p 0); 12M sig (prob 0.00203, p 0); 12M sig (prob 0.00166, p 0)
- FGF1_FGFR1: 3M sig (prob 0.00564, p 0); 3M sig (prob 0.00531, p 0); 3M sig (prob 0.00419, p 0); 3M sig (prob 0.00219, p 0.01); 3M sig (prob 0.00215, p 0)
- FGF1_FGFR2: 3M sig (prob 0.036, p 0); 3M sig (prob 0.0359, p 0); 3M sig (prob 0.0197, p 0); 3M sig (prob 0.0193, p 0); 3M sig (prob 0.0191, p 0); 3M sig (prob 0.00879, p 0.01); 3M sig (prob 0.00829, p 0.01)
- FN1_ITGAV_ITGB1: 3M sig (prob 0.0216, p 0); 3M sig (prob 0.016, p 0); 3M ns (prob 0.00816, p 0.42); 3M ns (prob 0.00627, p 0.43); 3M ns (prob 0.00602, p 0.47); 3M ns (prob 0.00581, p 0.45); 3M ns (prob 0.00536, p 0.43); 3M ns (prob 0.00463, p 0.52); 3M ns (prob 0.0044, p 0.44); 3M ns (prob 0.00429, p 0.55); 3M ns (prob 0.00395, p 0.51); 3M ns (prob 0.00326, p 0.45); 3M ns (prob 0.00324, p 0.56); 3M ns (prob 0.0024, p 0.56); 12M sig (prob 0.0157, p 0); 12M sig (prob 0.0114, p 0); 12M ns (prob 0.00547, p 0.32); 12M ns (prob 0.00523, p 0.31); 12M ns (prob 0.00394, p 0.68); 12M ns (prob 0.00377, p 0.99); 12M ns (prob 0.00365, p 0.56); 12M ns (prob 0.00291, p 0.6); 12M ns (prob 0.00263, p 1); 12M ns (prob 0.00209, p 1); 12M ns (prob 0.00152, p 0.61); 12M ns (prob 0.0011, p 1)
- GAS6_AXL: 3M sig (prob 0.00191, p 0)
- GAS6_MERTK: 3M sig (prob 0.0176, p 0); 3M sig (prob 0.00451, p 0); 3M sig (prob 0.00223, p 0); 3M sig (prob 0.00143, p 0)
- IGF1_IGF1R: 3M sig (prob 0.0116, p 0); 3M sig (prob 0.00983, p 0); 3M sig (prob 0.00768, p 0); 3M sig (prob 0.00697, p 0); 3M sig (prob 0.00684, p 0); 3M sig (prob 0.00503, p 0); 3M sig (prob 0.00473, p 0)
- LAMA2_ITGA6_ITGB1: 3M sig (prob 0.0169, p 0); 3M sig (prob 0.00855, p 0); 3M sig (prob 0.00713, p 0); 3M sig (prob 0.00533, p 0); 3M sig (prob 0.00465, p 0); 12M sig (prob 0.00693, p 0); 12M sig (prob 0.00617, p 0); 12M sig (prob 0.0026, p 0); 12M sig (prob 0.00205, p 0)
- LAMA4_ITGA6_ITGB1: 3M sig (prob 0.00258, p 0); 12M sig (prob 0.00273, p 0)
- LAMA5_SV2A: 3M sig (prob 0.0052, p 0)
- LAMB1_ITGA6_ITGB1: 3M sig (prob 0.00405, p 0); 12M sig (prob 0.00315, p 0)
- MPZ_MPZ: 3M sig (prob 0.00466, p 0)
- MPZ_MPZL1: 3M sig (prob 0.00669, p 0); 3M sig (prob 0.00618, p 0); 3M sig (prob 0.00466, p 0); 3M sig (prob 0.00333, p 0); 3M sig (prob 0.00256, p 0)
- NCAM1_NCAM1: 3M sig (prob 0.0625, p 0); 12M sig (prob 0.0127, p 0); 24M sig (prob 0.013, p 0)
- NRG1_ERBB4: 3M sig (prob 0.113, p 0); 3M sig (prob 0.0769, p 0); 3M sig (prob 0.075, p 0); 3M sig (prob 0.057, p 0); 3M sig (prob 0.0513, p 0); 3M sig (prob 0.0316, p 0); 3M sig (prob 0.0255, p 0); 24M sig (prob 0.043, p 0); 24M sig (prob 0.0285, p 0); 24M sig (prob 0.0273, p 0); 24M sig (prob 0.0265, p 0); 24M sig (prob 0.0161, p 0); 24M sig (prob 0.00981, p 0); 24M sig (prob 0.00873, p 0)
- NRG3_ERBB4: 3M sig (prob 0.275, p 0); 3M sig (prob 0.199, p 0); 3M sig (prob 0.195, p 0); 3M sig (prob 0.153, p 0); 3M sig (prob 0.139, p 0); 3M sig (prob 0.0888, p 0); 3M sig (prob 0.0724, p 0); 12M sig (prob 0.236, p 0); 12M sig (prob 0.18, p 0); 12M sig (prob 0.159, p 0); 12M sig (prob 0.141, p 0); 12M sig (prob 0.131, p 0); 12M sig (prob 0.0847, p 0); 12M sig (prob 0.0569, p 0); 24M sig (prob 0.283, p 0); 24M sig (prob 0.205, p 0); 24M sig (prob 0.198, p 0); 24M sig (prob 0.193, p 0); 24M sig (prob 0.126, p 0); 24M sig (prob 0.0802, p 0); 24M sig (prob 0.0719, p 0)
- NRXN1_NLGN1: 3M sig (prob 0.086, p 0); 3M sig (prob 0.0858, p 0); 3M sig (prob 0.0432, p 0); 3M sig (prob 0.0431, p 0); 3M sig (prob 0.0145, p 0); 3M sig (prob 0.007, p 0); 12M sig (prob 0.0761, p 0); 12M sig (prob 0.065, p 0); 12M sig (prob 0.00755, p 0); 24M sig (prob 0.0795, p 0); 24M sig (prob 0.0696, p 0); 24M sig (prob 0.0592, p 0); 24M sig (prob 0.0518, p 0)
- NRXN1_NLGN2: 3M sig (prob 0.0256, p 0); 3M sig (prob 0.0255, p 0); 3M ns (prob 0.00409, p 0.05); 12M sig (prob 0.0106, p 0); 12M sig (prob 0.00897, p 0); 12M sig (prob 0.00099, p 0)
- NRXN3_NLGN1: 3M sig (prob 0.0259, p 0); 3M sig (prob 0.013, p 0); 3M sig (prob 0.0126, p 0); 3M sig (prob 0.0124, p 0); 3M sig (prob 0.00913, p 0); 3M sig (prob 0.00627, p 0); 3M sig (prob 0.006, p 0); 3M sig (prob 0.0044, p 0.01); 12M sig (prob 0.0174, p 0); 12M sig (prob 0.00892, p 0); 12M sig (prob 0.00687, p 0); 24M sig (prob 0.0567, p 0); 24M sig (prob 0.042, p 0); 24M sig (prob 0.0131, p 0); 24M sig (prob 0.0121, p 0); 24M sig (prob 0.011, p 0); 24M sig (prob 0.0103, p 0); 24M sig (prob 0.00999, p 0); 24M sig (prob 0.00956, p 0); 24M sig (prob 0.00888, p 0); 24M sig (prob 0.00808, p 0); 24M sig (prob 0.00751, p 0); 24M sig (prob 0.0073, p 0)
- NTF3_NTRK2: 3M sig (prob 0.0366, p 0); 3M sig (prob 0.00291, p 0)
- NTF3_NTRK3: 3M sig (prob 0.0233, p 0)
- PROS1_AXL: 3M sig (prob 0.000621, p 0)
- PTN_ALK: 3M sig (prob 0.317, p 0); 3M sig (prob 0.274, p 0); 3M sig (prob 0.247, p 0); 3M sig (prob 0.239, p 0); 3M sig (prob 0.237, p 0); 3M sig (prob 0.233, p 0); 3M sig (prob 0.231, p 0); 12M sig (prob 0.21, p 0); 12M sig (prob 0.172, p 0); 12M sig (prob 0.17, p 0); 12M sig (prob 0.165, p 0); 12M sig (prob 0.162, p 0); 12M sig (prob 0.153, p 0); 12M sig (prob 0.147, p 0); 24M sig (prob 0.227, p 0); 24M sig (prob 0.168, p 0); 24M sig (prob 0.113, p 0); 24M sig (prob 0.0918, p 0); 24M sig (prob 0.0567, p 0); 24M sig (prob 0.0549, p 0); 24M sig (prob 0.0539, p 0)
- PTN_PTPRZ1: 3M sig (prob 0.217, p 0); 3M sig (prob 0.183, p 0); 3M sig (prob 0.163, p 0); 3M sig (prob 0.157, p 0); 3M sig (prob 0.156, p 0); 3M sig (prob 0.153, p 0); 3M sig (prob 0.152, p 0); 3M sig (prob 0.0756, p 0); 3M sig (prob 0.0622, p 0); 3M sig (prob 0.0547, p 0); 3M sig (prob 0.0524, p 0); 3M sig (prob 0.0519, p 0); 3M sig (prob 0.0508, p 0); 3M sig (prob 0.0502, p 0); 12M sig (prob 0.275, p 0); 12M sig (prob 0.228, p 0); 12M sig (prob 0.227, p 0); 12M sig (prob 0.221, p 0); 12M sig (prob 0.217, p 0); 12M sig (prob 0.205, p 0); 12M sig (prob 0.198, p 0); 12M sig (prob 0.133, p 0); 12M sig (prob 0.107, p 0); 12M sig (prob 0.106, p 0); 12M sig (prob 0.102, p 0); 12M sig (prob 0.1, p 0); 12M sig (prob 0.0939, p 0); 12M sig (prob 0.0903, p 0); 24M sig (prob 0.274, p 0); 24M sig (prob 0.206, p 0); 24M sig (prob 0.187, p 0); 24M sig (prob 0.141, p 0); 24M sig (prob 0.136, p 0); 24M sig (prob 0.115, p 0); 24M sig (prob 0.0908, p 0); 24M sig (prob 0.0731, p 0); 24M sig (prob 0.0715, p 0); 24M sig (prob 0.0692, p 0); 24M sig (prob 0.0679, p 0); 24M sig (prob 0.0448, p 0.02); 24M ns (prob 0.0434, p 0.05); 24M sig (prob 0.0426, p 0.04); 24M ns (prob 0.0354, p 0.06); 24M ns (prob 0.0246, p 0.11); 24M ns (prob 0.0157, p 0.17); 24M ns (prob 0.0124, p 0.17); 24M ns (prob 0.00744, p 0.19); 24M ns (prob 0.00719, p 0.2); 24M ns (prob 0.00705, p 0.19)
- PTN_SDC2: 3M sig (prob 0.0466, p 0); 3M ns (prob 0.0381, p 0.05); 3M ns (prob 0.0334, p 0.12); 3M ns (prob 0.032, p 0.17); 3M ns (prob 0.0316, p 0.15); 3M ns (prob 0.031, p 0.17); 3M ns (prob 0.0306, p 0.18); 3M sig (prob 0.0229, p 0.01); 3M ns (prob 0.0187, p 0.06); 3M ns (prob 0.0163, p 0.06); 3M ns (prob 0.0156, p 0.07); 3M ns (prob 0.0155, p 0.08); 3M ns (prob 0.0151, p 0.08); 3M ns (prob 0.0149, p 0.08); 12M sig (prob 0.0635, p 0); 12M sig (prob 0.0528, p 0); 12M sig (prob 0.0501, p 0); 12M sig (prob 0.0496, p 0); 12M sig (prob 0.048, p 0); 12M sig (prob 0.047, p 0); 12M sig (prob 0.0438, p 0); 12M sig (prob 0.042, p 0); 12M sig (prob 0.0416, p 0.03); 12M sig (prob 0.0412, p 0.03); 12M sig (prob 0.0399, p 0.03); 12M sig (prob 0.039, p 0.03); 12M ns (prob 0.0364, p 0.07); 12M ns (prob 0.0348, p 0.11); 12M ns (prob 0.0152, p 1); 12M ns (prob 0.0119, p 1); 12M ns (prob 0.0118, p 1); 12M ns (prob 0.0114, p 1); 12M ns (prob 0.0111, p 1); 12M ns (prob 0.0103, p 1); 12M ns (prob 0.0099, p 1); 24M sig (prob 0.0698, p 0); 24M sig (prob 0.049, p 0); 24M ns (prob 0.0374, p 0.07); 24M ns (prob 0.0316, p 0.25); 24M ns (prob 0.026, p 0.38); 24M ns (prob 0.0251, p 0.39); 24M ns (prob 0.0166, p 0.5); 24M ns (prob 0.0166, p 0.6); 24M ns (prob 0.0151, p 0.64); 24M ns (prob 0.0146, p 0.64); 24M ns (prob 0.0143, p 0.62); 24M ns (prob 0.0132, p 0.61); 24M ns (prob 0.0115, p 0.5); 24M ns (prob 0.00787, p 0.67); 24M ns (prob 0.0076, p 0.67); 24M ns (prob 0.00745, p 0.67); 24M ns (prob 0.0073, p 0.51); 24M ns (prob 0.00578, p 0.51); 24M ns (prob 0.00345, p 0.51); 24M ns (prob 0.00333, p 0.51); 24M ns (prob 0.00326, p 0.51)
- PTPRC_MRC1: 3M sig (prob 0.0643, p 0); 12M sig (prob 0.0883, p 0); 24M sig (prob 0.0447, p 0)
- SEMA3A_NRP1_PLXNA1: 3M sig (prob 0.0181, p 0); 3M sig (prob 0.00968, p 0); 3M sig (prob 0.0057, p 0); 3M sig (prob 0.00178, p 0)
- SEMA3A_NRP1_PLXNA2: 3M sig (prob 0.023, p 0); 3M sig (prob 0.021, p 0); 3M sig (prob 0.0123, p 0); 3M sig (prob 0.0113, p 0); 3M sig (prob 0.00727, p 0); 3M sig (prob 0.00664, p 0); 3M sig (prob 0.00228, p 0.03); 3M sig (prob 0.00208, p 0.04); 12M sig (prob 0.0141, p 0); 12M sig (prob 0.0126, p 0); 12M sig (prob 0.00815, p 0); 24M sig (prob 0.023, p 0); 24M sig (prob 0.0145, p 0); 24M sig (prob 0.00586, p 0)
- SPP1_ITGA9_ITGB1: 3M sig (prob 0.0222, p 0); 3M sig (prob 0.00503, p 0); 3M sig (prob 0.00485, p 0); 3M sig (prob 0.00339, p 0); 12M sig (prob 0.00708, p 0); 24M sig (prob 0.0041, p 0)
- SPP1_ITGAV_ITGB1: 3M sig (prob 0.0206, p 0); 3M sig (prob 0.0152, p 0); 3M sig (prob 0.00467, p 0.03); 3M sig (prob 0.0045, p 0.03); 3M sig (prob 0.00345, p 0); 3M sig (prob 0.00332, p 0); 3M ns (prob 0.00314, p 0.17); 3M ns (prob 0.00232, p 0.18); 12M sig (prob 0.0064, p 0); 12M sig (prob 0.00461, p 0.01)
- SPP1_ITGAV_ITGB5: 3M sig (prob 0.0232, p 0); 3M sig (prob 0.00528, p 0); 3M sig (prob 0.00509, p 0); 3M sig (prob 0.00355, p 0.03); 12M sig (prob 0.00774, p 0)
- TGFB2_TGFBR1_TGFBR2: 3M sig (prob 0.0235, p 0); 3M sig (prob 0.0135, p 0); 3M sig (prob 0.0134, p 0); 3M sig (prob 0.0087, p 0); 3M sig (prob 0.00494, p 0.04); 3M sig (prob 0.00493, p 0); 12M sig (prob 0.0193, p 0); 12M sig (prob 0.0105, p 0.01); 12M sig (prob 0.00712, p 0); 12M sig (prob 0.00394, p 0); 12M sig (prob 0.00232, p 0.02); 12M ns (prob 0.00211, p 0.35); 12M ns (prob 0.00144, p 0.33); 12M ns (prob 0.00124, p 0.34); 12M ns (prob 0.000846, p 0.31); 24M sig (prob 0.00868, p 0); 24M sig (prob 0.00301, p 0.01); 24M sig (prob 0.00258, p 0.01); 24M ns (prob 0.00162, p 0.08)
- VEGFA_VEGFR1: 3M sig (prob 0.0037, p 0); 3M sig (prob 0.00231, p 0); 3M sig (prob 0.00148, p 0)
