# iCluster 1 signature — external validation in MFS/UPS

This folder holds the code for the external validation of the iCluster 1 gene signature in myxofibrosarcoma (MFS) and undifferentiated pleomorphic sarcoma (UPS). The signature is tested in two public cohorts:

- **Mitra et al., *J Pathol* 2024** ([doi:10.1002/path.6347](https://doi.org/10.1002/path.6347)): n = 136 primary tumors, endpoint metastasis-free survival (MetFS).
- **TCGA SARC PanCancer Atlas**: n = 75 MFS/UPS, endpoints overall survival (OS) and disease-specific survival (DSS).

**Input data are not distributed with this repository.** Everything used is public. `scripts/01_download_data.R` downloads it, except for one file that has to be downloaded by hand (see [Data](#data)).

## Quick start

These steps use only R, so they work on Windows, macOS and Linux. All paths are relative to this `external_validation/` folder.

```r
# With external_validation/ as the working directory (in RStudio: Session > Set Working Directory)
source("scripts/00_install_packages.R")   # once
source("scripts/01_download_data.R")      # ~2.4 GB; then add Table S1 by hand (data/README.md)
source("run_all.R")                       # prepare data -> figures -> HR summary
```

The same steps from a terminal:

```bash
cd external_validation
Rscript scripts/00_install_packages.R
Rscript scripts/01_download_data.R
Rscript run_all.R
```

## Folder layout

```
├── R/                         shared functions (sourced by every script)
│   ├── setup.R                repo-relative paths, gene lists, package loading
│   ├── cohorts.R              cohort loading/filtering, GSVA scoring, median split
│   ├── plots.R                themes, CONSORT, lollipop, Kaplan–Meier, heatmaps, saving
│   └── stats.R                Cox + log-rank summary table
├── scripts/                   run in numeric order
│   ├── 00_install_packages.R  CRAN + Bioconductor dependencies
│   ├── 01_download_data.R     TCGA (cBioPortal) + E-MTAB-13947 (ArrayExpress)
│   ├── 02_prepare_data.R      Table S1 → CSV; RSEM isoforms → gene TPM; TCGA → RDS
│   ├── 03_mitra_2024.R        Mitra et al. (J Pathol 2024): Fig. 1D, Supp. Fig. 2B–D, MetFS HR
│   ├── 04_tcga_sarc.R         TCGA: Fig. 1E, Supp. Fig. 2E–G, OS/DSS HR
│   └── 05_survival_summary.R  HR table for both cohorts (xlsx)
├── run_all.R                  runs 02–05 (add `download` to also run 01)
├── data/                      inputs, not in git (see data/README.md)
├── figures/                   one PDF per manuscript panel
└── results/                   HR tables (CSV, xlsx)
```

## Figures → code

| Figure panel | Content | Script | Output (`figures/`) |
|---|---|---|---|
| Fig. 1D | Mitra et al., *J Pathol* 2024: MetFS Kaplan–Meier (HR 1.86, p = 0.034) | `03_mitra_2024.R` | [`Fig1D_Mitra2024_MetFS_KM_n136.pdf`](figures/Fig1D_Mitra2024_MetFS_KM_n136.pdf) |
| Fig. 1E | TCGA MFS/UPS: OS Kaplan–Meier (HR 2.29, p = 0.039) | `04_tcga_sarc.R` | [`Fig1E_TCGA_MFS_UPS_OS_KM_n75.pdf`](figures/Fig1E_TCGA_MFS_UPS_OS_KM_n75.pdf) |
| Supp. Fig. 2B | CONSORT diagram for the Mitra et al., *J Pathol* 2024 cohort (162 → 145 → 136) | `03_mitra_2024.R` | [`SuppFig2B_Mitra2024_CONSORT.pdf`](figures/SuppFig2B_Mitra2024_CONSORT.pdf) |
| Supp. Fig. 2C | Mitra et al., *J Pathol* 2024: ranked GSVA scores | `03_mitra_2024.R` | [`SuppFig2C_Mitra2024_GSVA_ranked_n136.pdf`](figures/SuppFig2C_Mitra2024_GSVA_ranked_n136.pdf) |
| Supp. Fig. 2D | Mitra et al., *J Pathol* 2024: signature heatmap | `03_mitra_2024.R` | [`SuppFig2D_Mitra2024_signature_heatmap_n136.pdf`](figures/SuppFig2D_Mitra2024_signature_heatmap_n136.pdf) |
| Supp. Fig. 2E | TCGA MFS/UPS: ranked GSVA scores | `04_tcga_sarc.R` | [`SuppFig2E_TCGA_MFS_UPS_GSVA_ranked_n75.pdf`](figures/SuppFig2E_TCGA_MFS_UPS_GSVA_ranked_n75.pdf) |
| Supp. Fig. 2F | TCGA MFS/UPS: signature heatmap with RPPA phospho-S6 | `04_tcga_sarc.R` | [`SuppFig2F_TCGA_MFS_UPS_signature_heatmap_n75.pdf`](figures/SuppFig2F_TCGA_MFS_UPS_signature_heatmap_n75.pdf) |
| Supp. Fig. 2G | TCGA MFS/UPS: DSS Kaplan–Meier (p = 0.091) | `04_tcga_sarc.R` | [`SuppFig2G_TCGA_MFS_UPS_DSS_KM_n75.pdf`](figures/SuppFig2G_TCGA_MFS_UPS_DSS_KM_n75.pdf) |

Each script writes a PDF and a 300 dpi PNG of every panel. Only the PDFs are committed; the PNGs are local previews and are excluded by `.gitignore`.

The R figures label the two sample groups "iCluster 1-high" and "iCluster 1-low". Supplementary Fig. 2 in the manuscript calls the same groups "iCluster 1-like" and "non-iCluster 1". The panels were assembled into the final figures outside R.

## Methods summary

- **Signature:** iCluster 1 genes: *TNC, TMEM158, STEAP1, SERPINE1, ITGA10, LRRC15*. They are scored per sample with GSVA (`kcdf = "Gaussian"`) on log2(TPM + 1) for Mitra et al. and log2(RSEM + 1) for TCGA.
- **Grouping:** samples at or above the cohort median GSVA score are iCluster 1-high (or "-like"); the rest are iCluster 1-low.
- **Survival testing:** log-rank test, then a Cox model with iCluster 1-low as the reference group.
- **Heatmaps:** row z-scores for the iCluster 1 and iCluster 2 genes (*CLEC3B, TNXB, MFAP5, GPR64, ABCA8*). In TCGA, RPPA S6, S6 pS235/S236 and S6 pS240/S244 are added below. *GPR64* is annotated as *ADGRG2* in current Ensembl and is renamed back to *GPR64*.
- **Mitra et al., *J Pathol* 2024 cohort filtering:** start with 162 MFS/UPS. Remove 17 local recurrences (LR1/LR2/LR3), leaving 145 primaries. Remove 3 with synchronous metastases (MFS time = 0) and 6 with no MFS follow-up, leaving 136. GSVA is computed on all 162 samples; the median is taken within the 136.
- **TCGA cohort:** ONCOTREE codes `MFH` (UPS) and `MFS`, giving n = 75.

### Hazard ratios

| Cohort | Endpoint | N | Events | HR (95% CI) | Log-rank p |
|---|---|---|---|---|---|
| Mitra et al., *J Pathol* 2024 | MetFS | 136 | 48 | 1.86 (1.04–3.32) | 0.034 |
| TCGA | OS | 75 | 25 | 2.29 (1.02–5.13) | 0.039 |
| TCGA | DSS | 75 | 19 | 2.16 (0.87–5.41) | 0.091 |

HR compares iCluster 1-high with iCluster 1-low. Full output is in `results/survival_HR_pvalues_Mitra_TCGA.xlsx`.

## Data

Where each file goes, plus the manual download steps, are in **[data/README.md](data/README.md)**.

### TCGA SARC PanCancer Atlas 2018 (cBioPortal)
- **Study:** https://www.cbioportal.org/study/summary?id=sarc_tcga_pan_can_atlas_2018
- **Files used:** `data_mrna_seq_v2_rsem.txt`, `data_rppa.txt`, `data_clinical_patient.txt`, `data_clinical_sample.txt`. They are downloaded automatically from the [cBioPortal datahub](https://github.com/cBioPortal/datahub/tree/master/public/sarc_tcga_pan_can_atlas_2018).
- **Citations:** Hoadley KA et al., *Cell* 2018;173:291–304; Cerami E et al., *Cancer Discov* 2012;2:401–404; Gao J et al., *Sci Signal* 2013;6:pl1.

### Mitra et al., *J Pathol* 2024: MFS/UPS cohort
> Mitra S, Farswan A, Piccinelli P, Sydow S, Hesla A, Tsagkozis P, Vult von Steyern F, Almqvist M, Eriksson M, Magnusson L, Nilsson J, Pillay N, Mertens F. Transcriptomic profiles of myxofibrosarcoma and undifferentiated pleomorphic sarcoma correlate with clinical and genomic features. *J Pathol* 2024;264:293–304. https://doi.org/10.1002/path.6347

The authors' data availability statement:

> Global gene expression values from 162 cases of MFS and undifferentiated pleomorphic sarcoma UPS, as well as from 21 control tumours, and copy number data on 145 cases of MFS and UPS have been uploaded to Annotare (https://www.ebi.ac.uk/fg/annotare/) with ArrayExpress accession numbers E-MTAB-13947 and E-MTAB-13948 respectively. The analysis of SNVs was performed by a commercial company (Eurofins, https://go.eurofinsgenomics.eu/). The unfiltered reports from the company could be obtained from the corresponding author upon reasonable request.

| Dataset | Source | Used here | How to get it |
|---|---|---|---|
| Transcriptomic analysis of 162 myxofibrosarcomas and undifferentiated pleomorphic sarcomas | [E-MTAB-13947](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13947) | **Yes.** The 162 tumor RSEM isoform files; the 21 controls are not used. | Automatic |
| Segment files from SNP array analysis of 145 myxofibrosarcoma and undifferentiated pleomorphic sarcoma | [E-MTAB-13948](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13948) | No. Not needed for these figures. | — |
| Table S1: clinical annotations and follow-up | Supporting Information, `path6347-sup-0002-tables1.xlsx` | **Yes.** Gives diagnosis, sample type and MFS time/event. | **Manual** (the publisher blocks scripted downloads). See [data/README.md](data/README.md). |

## Software

These figures were made with R 4.5.0 and the following packages: data.table 1.18.4, dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3, patchwork 1.3.2, cowplot 1.1.3, survival 3.8.9, survminer 0.5.0, readxl 1.4.5, writexl 1.5.4. From Bioconductor: GSVA 2.2.0, ComplexHeatmap 2.24.0, circlize 0.4.16, org.Hs.eg.db 3.21.0, AnnotationDbi 1.70.0.
