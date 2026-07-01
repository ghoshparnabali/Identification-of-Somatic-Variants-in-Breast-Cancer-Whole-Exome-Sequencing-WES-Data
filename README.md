# Identification of Somatic Variants in Breast Cancer Whole-Exome Sequencing (WES) Data
This repository contains an end-to-end **tumour-only somatic variant-calling pipeline** built in Python (Jupyter / Google Colab) with Visual Studio Code to identify candidate somatic mutations in **triple-negative breast cancer (TNBC)** from whole-exome sequencing data. Using tumour exomes from the 2026 proteogenomic study [PRJNA1422845](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1422845) (SRA study [SRP676624](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP676624&o=acc_s%3Aa)), conducted by Yonsei Cancer Center, Yonsei University College of Medicine, the pipeline follows the GATK Best Practices workflow for somatic short-variant discovery — BWA-MEM alignment, duplicate marking, base-quality recalibration, and Mutect2 tumour-only calling with a germline resource, panel of normals, and orientation-bias modelling — followed by snpEff functional annotation and prioritisation against a hereditary breast-cancer gene panel. The reference genome is subset to chromosomes 13 and 17 to focus the analysis on **BRCA2, BRCA1, and TP53** while remaining computationally feasible on a Colab CPU runtime.

---

## Key Technical Implementations
- **Tumour-only somatic calling with GATK Mutect2 —** The dataset provides tumour exomes without matched normal samples, so a conventional tumour–normal subtraction is not possible. Instead, a population germline resource (gnomAD allele frequencies) and a Panel of Normals stand in for the matched normal: they let Mutect2 estimate the probability that each variant is germline or a recurrent technical artifact rather than a true somatic event. Because no matched normal is available, the output is reported as *candidate* somatic variants — rare germline variants cannot be formally excluded.
- **Two-chromosome reference subset (chr13 + chr17) —** Rather than aligning against the whole genome, the reference is restricted to chromosome 13 (which carries *BRCA2*) and chromosome 17 (which carries *BRCA1* and *TP53*) — the core genes of hereditary and sporadic breast cancer. This is a deliberate design choice that keeps alignment and variant calling tractable on Colab's limited CPU while covering the genes of interest. The trade-off is that some reads originating elsewhere in the genome are force-aligned to these two chromosomes, producing clustered false calls, which are addressed downstream by Mutect2's filters and by prioritisation to known cancer genes.
- **Orientation-bias and contamination modelling —** Somatic mutations often occur at low allele fractions, where sequencing artifacts can masquerade as real variants. The pipeline runs `LearnReadOrientationModel` to model strand-orientation artifacts (e.g. oxidative/FFPE damage) and `GetPileupSummaries` / `CalculateContamination` to estimate cross-sample contamination, both of which are supplied to `FilterMutectCalls` so that low-confidence calls are flagged rather than reported as findings.
- **Base Quality Score Recalibration (BQSR) —** Sequencers introduce systematic, machine-specific errors into base-quality scores. Using known-variant sites (dbSNP, Mills indels), BQSR recalibrates these scores empirically so that Mutect2's probabilistic model operates on accurate confidence values, reducing false-positive calls.
- **Rarity-based germline exclusion and local annotation with snpEff —** Because a matched normal is unavailable, population allele frequency does double duty as a germline filter: variants that are common in gnomAD are treated as likely germline, while rare or absent variants are retained as candidate somatic. gnomAD allele frequencies are transferred onto the called variants with `bcftools annotate`, and snpEff provides gene-level functional consequence, protein change, and predicted impact — annotating locally with no dependence on external servers.

## Key Findings
- **Candidate somatic loss-of-function mutations in TP53 were identified in 2 of the 3 tumours.** Both are protein-truncating (a nonsense mutation and a frameshift), classified as HIGH impact, and absent from gnomAD — a pattern consistent with somatic driver mutations. TP53 is the most frequently mutated gene in triple-negative breast cancer and is characteristically hit by truncating loss-of-function mutations, so these findings reflect the expected TNBC driver biology and validate the pipeline end to end.

### Cancer-panel candidate variants
| Sample | Gene | Position (GRCh38) | Change | Consequence | Impact | gnomAD |
|--------|------|-------------------|--------|-------------|--------|--------|
| Sample_1 | **TP53** | chr17:7670664 | p.Glu310* | stop_gained | HIGH | absent |
| Sample_2 | **TP53** | chr17:7676045 | p.Gly69fs | frameshift_variant | HIGH | absent |

- **Tumour heterogeneity is evident across the cohort.** While two tumours carried candidate somatic *TP53* mutations, the third had no candidate mutation in the surveyed cancer genes — a legitimate and unsurprising result given that TP53, though dominant in TNBC, is not universal.

### Per-sample candidate summary
| Sample | SRA run | Candidate variants | In cancer panel |
|--------|---------|--------------------|-----------------|
| Sample_1 | SRR37211323 | 174 | 1 (*TP53*) |
| Sample_2 | SRR37211334 | 169 | 1 (*TP53*) |
| Sample_3 | SRR37211335 | 73 | 0 |

- **Interpretation and scope.** The candidate variants are reported as *candidate somatic* because the tumour-only design cannot exclude rare germline variants without a matched normal. The larger per-sample candidate counts are dominated by novel/uncharacterised variants and include alignment artifacts inherent to the two-chromosome subset reference; the cancer-panel hits above are the biologically interpretable findings. Full quality-control reports (pre- and post-trim) are provided in `multiqc_report/`, and the complete candidate list and summary report are in `results/`.

---

## Repository Structure

```
Identification_of_Somatic_Variants_in_Breast_Cancer_Whole_Exome_Sequencing_WES_Data/
├── WES_Variant_Analysis/
│   ├── code/
│   │   └── WES_Variant_Analysis.ipynb
│   ├── multiqc_report/
│   │   ├── multiqc_report_pre_trim.html
│   │   └── multiqc_report_post_trim.html
│   └── results/
│       ├── analysis_report.md
│       └── candidate_somatic_variants.csv
├── LICENSE
└── README.md
```

---

## Data Acquisition

The analysis uses three tumour whole-exome samples from NCBI BioProject [PRJNA1422845](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1422845) (SRA study [SRP676624](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP676624&o=acc_s%3Aa)), a 2026 proteogenomic study of triple-negative breast cancer. The data are openly available (public consent). The three smallest tumour runs were selected:

| Sample | SRA run accession |
|--------|-------------------|
| Sample_1 | SRR37211323 |
| Sample_2 | SRR37211334 |
| Sample_3 | SRR37211335 |

Raw sequencing data (FASTQ), aligned BAM files, the reference genome, and the GATK resource bundles are **not** included in this repository due to their size; they are downloaded and generated automatically when the notebook is run. The pipeline retrieves the raw reads directly from SRA using the accessions above.

---

## Dependencies

All analyses were performed on a **Google Colab** CPU runtime (Ubuntu), driven through VS Code. The pipeline uses the following tools, installed within the notebook:

| Tool | Version | Purpose |
|------|---------|---------|
| SRA-Toolkit | latest | Downloading raw reads from SRA |
| FastQC | latest | Per-sample read quality control |
| MultiQC | latest | Aggregated QC reporting |
| Cutadapt | latest | Adapter and quality trimming |
| BWA | latest | Read alignment (BWA-MEM) |
| SAMtools | 1.13 | BAM sorting, indexing, statistics |
| BCFtools | 1.13 | VCF filtering and annotation transfer |
| GATK | 4.5.0.0 | MarkDuplicates, BQSR, Mutect2, filtering |
| snpEff | v4_3 (GRCh38.86 database) | Functional variant annotation |
| OpenJDK | 17 | Java runtime for GATK and snpEff |

Reference and resource files (GATK Best Practices, hg38): the GRCh38 reference (chr13 + chr17 subset), af-only-gnomAD germline resource, 1000G Panel of Normals, and known-sites for BQSR (dbSNP, Mills indels) are downloaded and subset within the notebook.

A representative install block (as run at the top of the notebook):

```bash
apt-get -qq install -y sra-toolkit fastqc bwa samtools bcftools tabix openjdk-17-jdk-headless
pip -q install cutadapt multiqc
# GATK 4.5.0.0 and snpEff are downloaded via wget within the notebook
```

---

## How to Reproduce

1. Clone this repository:
   ```
   git clone https://github.com/ghoshparnabali/Identification-of-Somatic-Variants-in-Breast-Cancer-Whole-Exome-Sequencing-WES-Data.git
   ```
2. Open `WES_Variant_Analysis/code/WES_Variant_Analysis.ipynb` in Google Colab (or a Jupyter environment with an Ubuntu backend).
3. Run the notebook cells from top to bottom. The pipeline will, in order: install the required tools, build and index the chr13 + chr17 reference, download and subset the GATK resource bundles, retrieve the three tumour exomes from SRA, and run quality control, trimming, alignment, duplicate marking, BQSR, Mutect2 tumour-only calling, contamination and orientation-bias filtering, snpEff annotation, and variant prioritisation.
4. Outputs are written to the working directory: quality-control reports, the annotated and filtered VCFs, the prioritised candidate list (`candidate_somatic_variants.csv`), and a summary report (`analysis_report.md`).

> **Note:** The alignment step is the most time-consuming stage on a Colab CPU runtime. Because the runtime disk is ephemeral, backing up intermediate BAM files to Google Drive is recommended to guard against session disconnects.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
