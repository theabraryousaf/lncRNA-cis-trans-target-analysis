#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
source config/config.env
mkdir -p data ref "$OUTDIR"

echo "Step 1: Save software versions"
{
  echo "Python:"
  python --version
  echo "R:"
  Rscript --version
  echo "BEDTools:"
  bedtools --version
} > "$OUTDIR/software_versions.txt"
cat "$OUTDIR/software_versions.txt"

echo "Step 2: Download GENCODE if needed"
bash scripts/01_download_gencode.sh

echo "Step 3: Validate input"
python scripts/02_validate_input.py \
  --input "$INPUT_CSV" \
  --sample-columns "$SAMPLE_COLUMNS"

echo "Step 4: Annotate genes and prepare BED files"
python scripts/03_annotate_and_prepare_bed.py \
  --input "$INPUT_CSV" \
  --gtf "$GTF_FILE" \
  --outdir "$OUTDIR" \
  --sample-columns "$SAMPLE_COLUMNS" \
  --lnc-biotypes "$LNC_BIOTYPES" \
  --pcg-biotypes "$PCG_BIOTYPES"

echo "Step 5: Run cis analysis"
bash scripts/04_run_cis_bedtools.sh

python scripts/05_parse_cis_targets.py \
  --raw "$OUTDIR/cis_raw_lncRNA_PCG_${CIS_WINDOW_BP}bp.tsv" \
  --outdir "$OUTDIR" \
  --window-bp "$CIS_WINDOW_BP"

echo "Step 6: Create log2 expression matrix"
python scripts/06_make_log2_expression_matrix.py \
  --input "$INPUT_CSV" \
  --sample-columns "$SAMPLE_COLUMNS" \
  --output "$OUTDIR/log2_normalized_expression_matrix.tsv"

echo "Step 7: Run trans analysis"
Rscript scripts/07_trans_targets_pearson.R \
  "$OUTDIR/log2_normalized_expression_matrix.tsv" \
  "$OUTDIR/significant_genes_annotated.tsv" \
  "$OUTDIR" \
  "$PCC_CUTOFF" \
  "$FDR_CUTOFF"

echo "Step 8: Combine cis and trans evidence"
python scripts/08_combine_cis_trans.py

echo "Step 9: Create summary report"
python scripts/09_make_summary_report.py

echo "Done. Main outputs are in: $OUTDIR"
echo "Important files:"
echo "  $OUTDIR/FINAL_lncRNA_cis_trans_priority.csv"
echo "  $OUTDIR/cis_targets_${CIS_WINDOW_BP}bp_all.csv"
echo "  $OUTDIR/trans_targets_Pearson_selected.csv"
echo "  $OUTDIR/all_cis_target_genes.csv"
echo "  $OUTDIR/all_trans_target_genes.csv"
echo "  $OUTDIR/FINAL_cis_trans_analysis_summary.txt"
