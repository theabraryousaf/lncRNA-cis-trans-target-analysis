# lncRNA cis/trans target analysis for macOS M-chip

A simple beginner-friendly pipeline to predict **candidate cis- and trans-target protein-coding genes of lncRNAs** from an RNA-seq normalized-count matrix.

This repository was designed for Apple Silicon Macs (`M1`, `M2`, `M3`, `M4`) and avoids FEELnc installation problems by using:

- **BEDTools** for cis-target prediction
- **R Pearson correlation** for trans-target prediction
- **GENCODE GTF annotation** for gene biotype, gene symbol and genomic coordinates

> Important: This pipeline predicts **candidate associations**, not experimentally validated regulatory targets.

---

## What this pipeline does

Given a CSV file containing significant genes from DESeq2 or another differential-expression analysis:

```text
gene_id,0h_1,0h_2,0h_3,16h_1,16h_2,16h_3,20h_1,20h_2,20h_3,24h_1,24h_2,24h_3
ENSG00000000457.14,646.58,685.77,783.36,...
```

The pipeline will:

1. Validate the input file.
2. Download or use a GENCODE GTF annotation.
3. Annotate each gene as `lncRNA`, `protein_coding` or other biotype.
4. Extract gene-level BED files for lncRNAs and protein-coding genes.
5. Predict **cis targets**: protein-coding genes within ±100 kb of each lncRNA.
6. Predict **trans targets**: protein-coding genes correlated with lncRNAs across samples.
7. Keep both positive and negative correlations using `abs(PCC)`.
8. Create final ranked candidate lncRNA tables.
9. Create final lists of all cis-target genes and all trans-target genes.

---

## Folder structure

```text
lncrna-cis-trans-macosm/
├── config/
│   └── config.env
├── data/
│   ├── README.md
│   └── example_input_template.csv
├── docs/
│   └── method_text.md
├── ref/
├── results/
├── scripts/
│   ├── 00_install_macos_m.sh
│   ├── 01_download_gencode.sh
│   ├── 02_validate_input.py
│   ├── 03_annotate_and_prepare_bed.py
│   ├── 04_run_cis_bedtools.sh
│   ├── 05_parse_cis_targets.py
│   ├── 06_make_log2_expression_matrix.py
│   ├── 07_trans_targets_pearson.R
│   ├── 08_combine_cis_trans.py
│   ├── 09_make_summary_report.py
│   └── run_all.sh
├── environment-macos-m.yml
├── .gitignore
└── README.md
```

---

## Step 1: Install tools on macOS M-chip

Open Terminal and run:

```bash
brew install micromamba git wget curl bedtools
```

Initialize micromamba for zsh if you have not already done it:

```bash
micromamba shell init -s zsh -r ~/micromamba
source ~/.zshrc
```

Create the environment:

```bash
cd ~/Documents/lncrna-cis-trans-macosm
micromamba create -y -n lncrna_cistrans -f environment-macos-m.yml
micromamba activate lncrna_cistrans
```

Check versions:

```bash
python --version
Rscript --version
bedtools --version
```

---

## Step 2: Prepare your input file

Put your CSV file inside the `data/` folder.

Example:

```bash
cp ~/Downloads/input_p0.05_lncrna_pcg.csv data/input_p0.05_lncrna_pcg.csv
```

Your file must contain:

- one `gene_id` column with Ensembl gene IDs
- normalized counts for each sample
- at least 3 sample columns

Example:

```text
gene_id,0h_1,0h_2,0h_3,16h_1,16h_2,16h_3,20h_1,20h_2,20h_3,24h_1,24h_2,24h_3
ENSG00000000457.14,646.58,685.77,783.36,1042.08,982.60,1038.50,1092.22,1107.30,1174.00,1553.01,1132.40,1330.04
```

The normalized counts should **not** be raw counts. DESeq2 normalized counts, TPM, CPM or another normalized expression scale can be used. The pipeline transforms values as `log2(count + 1)` before correlation.

---

## Step 3: Edit the configuration file

Open:

```bash
nano config/config.env
```

Important settings:

```bash
INPUT_CSV="data/input_p0.05_lncrna_pcg.csv"
SAMPLE_COLUMNS="0h_1,0h_2,0h_3,16h_1,16h_2,16h_3,20h_1,20h_2,20h_3,24h_1,24h_2,24h_3"
CIS_WINDOW_BP="100000"
PCC_CUTOFF="0.85"
FDR_CUTOFF="0.01"
```

For trans analysis:

- `PCC_CUTOFF=0.85` keeps strong positive and negative correlations.
- `FDR_CUTOFF=0.01` keeps statistically stronger pairs.
- The original `PCC` sign is retained to classify targets as positive or negative.

---

## Step 4: Run the full analysis

```bash
cd ~/Documents/lncrna-cis-trans-macosm
micromamba activate lncrna_cistrans
bash scripts/run_all.sh
```

---

## Main output files

### Final ranked lncRNAs

```text
results/FINAL_lncRNA_cis_trans_priority.csv
results/FINAL_lncRNA_cis_trans_priority.tsv
```

### All cis pairs

```text
results/cis_targets_100kb_all.csv
results/cis_targets_100kb_all.tsv
```

### All selected trans pairs

```text
results/trans_targets_Pearson_selected.csv
results/trans_targets_Pearson_selected.tsv
```

### All unique cis-target protein-coding genes

```text
results/all_cis_target_genes.csv
results/all_cis_target_genes.txt
```

### All unique trans-target protein-coding genes

```text
results/all_trans_target_genes.csv
results/all_trans_target_genes.txt
```

### Positive and negative trans-target genes

```text
results/all_trans_positive_target_genes.csv
results/all_trans_positive_target_genes.txt
results/all_trans_negative_target_genes.csv
results/all_trans_negative_target_genes.txt
```

### Summary report

```text
results/FINAL_cis_trans_analysis_summary.txt
```

---

## How to interpret PCC

PCC means **Pearson correlation coefficient**.

- `PCC > 0`: lncRNA and protein-coding gene increase/decrease together.
- `PCC < 0`: lncRNA and protein-coding gene move in opposite directions.
- `abs(PCC)` is used for selection so both positive and negative associations are retained.

Example:

```text
PCC = +0.91  positive candidate association
PCC = -0.91  negative candidate association
```

Use `PCC` for direction and `abs_PCC` for strength.

---

## Recommended manuscript wording

Candidate cis-target genes were defined as protein-coding genes located within 100 kb upstream or downstream of each lncRNA locus. Candidate trans-target genes were identified by Pearson correlation analysis between lncRNA and protein-coding gene expression profiles across all samples. Pairs with absolute Pearson correlation coefficient ≥0.85 and Benjamini–Hochberg FDR ≤0.01 were retained, and the original PCC sign was used to distinguish positive and negative associations.

---

## Notes for beginners

- Use the same GTF annotation that was used for your RNA-seq quantification if possible.
- Do not use mean normalized counts alone for trans analysis. You need sample-wise expression values.
- Do not call these validated targets. They are candidate cis/trans targets.
- Change `PCC_CUTOFF` and `FDR_CUTOFF` in `config/config.env` if you need looser or stricter filtering.
