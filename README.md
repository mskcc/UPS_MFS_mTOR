# Li, G, et. al., Molecular Cancer Therapeutics 2026

## Genomic clustering identifies an aggressive subset of myxofibrosarcoma and undifferentiated pleomorphic sarcoma sensitive to combined mTOR and MEK inhibition

Computational methods for MCT-26-0243-ETR

## Contents

| Folder | Description |
|---|---|
| [`external_validation/`](external_validation/) | Validation of the iCluster 1 signature in two public MFS/UPS cohorts (Mitra et al., *J Pathol* 2024; TCGA SARC), with the code for Fig. 1D–E and Supplementary Fig. 2B–G |

---

## External validation of the iCluster 1 signature

We tested whether the iCluster 1 gene signature predicts outcome in myxofibrosarcoma (MFS) and undifferentiated pleomorphic sarcoma (UPS) outside our cohort. It was scored in two public cohorts that were not used to derive it.

### Cohorts

| Cohort | Data | Patients | Endpoint |
|---|---|---|---|
| **Mitra et al., *J Pathol* 2024** ([doi:10.1002/path.6347](https://doi.org/10.1002/path.6347)) | RNA-seq: [E-MTAB-13947](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13947). Clinical data: Table S1 of the paper. | 136 primary MFS/UPS | Metastasis-free survival (MetFS) |
| **TCGA SARC PanCancer Atlas 2018** ([cBioPortal](https://www.cbioportal.org/study/summary?id=sarc_tcga_pan_can_atlas_2018)) | RNA-seq (RSEM), RPPA, clinical | 75 MFS/UPS (ONCOTREE `MFS`, `MFH`) | Overall survival (OS) and disease-specific survival (DSS) |

The Mitra et al., *J Pathol* 2024 cohort starts with 162 MFS/UPS cases. Removing 17 local recurrences leaves 145 primary tumors. Removing 3 with synchronous metastases and 6 without MFS follow-up leaves the 136 analyzed (Supp. Fig. 2B).

### Approach

1. **Signature score:** each tumor is scored for the six iCluster 1 genes (*TNC, TMEM158, STEAP1, SERPINE1, ITGA10, LRRC15*) with GSVA on log2-transformed expression.
2. **Grouping:** within each cohort, tumors at or above the median score are iCluster 1-high (or "iCluster 1-like"); the rest are iCluster 1-low.
3. **Survival testing:** Kaplan–Meier curves, the log-rank test, and a Cox model with iCluster 1-low as the reference group.
4. **Heatmaps:** heatmaps of the iCluster 1 and iCluster 2 marker genes (*CLEC3B, TNXB, MFAP5, GPR64, ABCA8*) confirm that the score separates the expected expression programs. In TCGA, RPPA phospho-S6 (S235/S236 and S240/S244), a readout of mTOR pathway activity, is shown alongside.

### Results

| Cohort | Endpoint | N | Events | HR (95% CI) | Log-rank p |
|---|---|---|---|---|---|
| Mitra et al., *J Pathol* 2024 | MetFS | 136 | 48 | 1.86 (1.04–3.32) | 0.034 |
| TCGA SARC | OS | 75 | 25 | 2.29 (1.02–5.13) | 0.039 |
| TCGA SARC | DSS | 75 | 19 | 2.16 (0.87–5.41) | 0.091 |

HR compares iCluster 1-high with iCluster 1-low. iCluster 1-high tumors had shorter metastasis-free survival in the Mitra et al. cohort and shorter overall survival in TCGA. DSS in TCGA trended the same way but did not reach significance.

### Figures

| Manuscript panel | Content | File |
|---|---|---|
| Fig. 1D | Mitra et al., *J Pathol* 2024: metastasis-free survival Kaplan–Meier (n = 136) | [`Fig1D_Mitra2024_MetFS_KM_n136.pdf`](external_validation/figures/Fig1D_Mitra2024_MetFS_KM_n136.pdf) |
| Fig. 1E | TCGA MFS/UPS: overall survival Kaplan–Meier (n = 75) | [`Fig1E_TCGA_MFS_UPS_OS_KM_n75.pdf`](external_validation/figures/Fig1E_TCGA_MFS_UPS_OS_KM_n75.pdf) |
| Supp. Fig. 2B | CONSORT diagram for the Mitra et al., *J Pathol* 2024 cohort | [`SuppFig2B_Mitra2024_CONSORT.pdf`](external_validation/figures/SuppFig2B_Mitra2024_CONSORT.pdf) |
| Supp. Fig. 2C | Mitra et al., *J Pathol* 2024: ranked iCluster 1 GSVA scores with median split | [`SuppFig2C_Mitra2024_GSVA_ranked_n136.pdf`](external_validation/figures/SuppFig2C_Mitra2024_GSVA_ranked_n136.pdf) |
| Supp. Fig. 2D | Mitra et al., *J Pathol* 2024: iCluster 1 / iCluster 2 gene heatmap | [`SuppFig2D_Mitra2024_signature_heatmap_n136.pdf`](external_validation/figures/SuppFig2D_Mitra2024_signature_heatmap_n136.pdf) |
| Supp. Fig. 2E | TCGA MFS/UPS: ranked iCluster 1 GSVA scores with median split | [`SuppFig2E_TCGA_MFS_UPS_GSVA_ranked_n75.pdf`](external_validation/figures/SuppFig2E_TCGA_MFS_UPS_GSVA_ranked_n75.pdf) |
| Supp. Fig. 2F | TCGA MFS/UPS: gene heatmap with RPPA phospho-S6 | [`SuppFig2F_TCGA_MFS_UPS_signature_heatmap_n75.pdf`](external_validation/figures/SuppFig2F_TCGA_MFS_UPS_signature_heatmap_n75.pdf) |
| Supp. Fig. 2G | TCGA MFS/UPS: disease-specific survival Kaplan–Meier | [`SuppFig2G_TCGA_MFS_UPS_DSS_KM_n75.pdf`](external_validation/figures/SuppFig2G_TCGA_MFS_UPS_DSS_KM_n75.pdf) |

Each panel is saved as its own PDF in `external_validation/figures/`. The panels were assembled into the final figures outside R.

### Reproducibility

The analysis is entirely in R and runs on Windows, macOS or Linux. The input data are public and are not stored in this repository. A download script fetches them, except for Table S1 of Mitra et al., *J Pathol* 2024, which has to be downloaded by hand ([instructions](external_validation/data/README.md)).

```bash
cd external_validation
Rscript scripts/00_install_packages.R   # CRAN + Bioconductor dependencies (once)
Rscript scripts/01_download_data.R      # TCGA + E-MTAB-13947 (~2.4 GB)
Rscript run_all.R                       # prepare data, make figures, write HR tables
```

[`external_validation/README.md`](external_validation/README.md) has the full methods, folder layout, data sources and package versions.
