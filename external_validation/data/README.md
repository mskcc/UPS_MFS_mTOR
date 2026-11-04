# Data

This repository does not include any input data. Everything used is public. `scripts/01_download_data.R` downloads all of it except one file, Table S1 of Mitra et al., *J Pathol* 2024, which has to be downloaded by hand.

```
data/
├── raw/                                   ← inputs (not in git)
│   ├── sarc_tcga_pan_can_atlas_2018/      automatic
│   │   ├── data_mrna_seq_v2_rsem.txt
│   │   ├── data_rppa.txt
│   │   ├── data_clinical_patient.txt
│   │   └── data_clinical_sample.txt
│   └── mitra_2024/
│       ├── path6347-sup-0002-tables1.xlsx MANUAL (see below)
│       └── E-MTAB-13947/                  automatic
│           ├── E-MTAB-13947.sdrf.txt
│           └── P*_fastq_Aligned.toTranscriptome.out.isoforms.results   (162 files)
└── processed/                             ← built by scripts/02_prepare_data.R (not in git)
    ├── mitra_2024_table_s1.csv
    ├── mitra_2024_tpm_matrix.rds
    └── tcga_sarc_data.rds
```

## Manual download: Mitra et al., *J Pathol* 2024, Table S1

The publisher (Wiley) blocks scripted downloads, so this one file has to be downloaded in a browser:

1. Open the article at https://doi.org/10.1002/path.6347.
2. Scroll to **Supporting Information**.
3. Download **`path6347-sup-0002-tables1.xlsx`** ("Table S1. General data on cohort of transcriptomically analysed myxofibrosarcomas and undifferentiated pleomorphic sarcomas").
4. Save it unchanged, with the same filename, at `data/raw/mitra_2024/path6347-sup-0002-tables1.xlsx`.

The analysis uses these Table S1 columns:

| Column | Used for |
|---|---|
| `Pat ID` | Patient ID |
| `Diagnosisa` | MFS or UPS |
| `Samplef` | Sample type: P = primary, LR1/LR2/LR3 = local recurrence |
| `MFS_Timee` | Metastasis-free survival time, in months |
| `Fud` | Follow-up text. A sample counts as an event if the text contains "Met". |

## Manual fallback for the automatic files

If `01_download_data.R` fails (for example behind a firewall), download the files in a browser and save them at the paths shown above.

- **TCGA SARC PanCancer Atlas 2018:** go to https://www.cbioportal.org/study/summary?id=sarc_tcga_pan_can_atlas_2018 and click the download icon. Unzip the archive and copy only the 4 files listed above.
- **E-MTAB-13947:** go to https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13947. Under "Files", download the 162 tumor files named `P<number>_fastq_Aligned.toTranscriptome.out.isoforms.results`. The `Control_*` files are not needed. They can also be downloaded straight from https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/947/E-MTAB-13947/Files/.
