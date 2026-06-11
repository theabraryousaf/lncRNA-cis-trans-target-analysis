#!/usr/bin/env bash
set -euo pipefail

source config/config.env
mkdir -p "$OUTDIR"

bedtools window \
  -w "$CIS_WINDOW_BP" \
  -a "$OUTDIR/significant_lncRNAs.bed" \
  -b "$OUTDIR/significant_protein_coding.bed" \
  > "$OUTDIR/cis_raw_lncRNA_PCG_${CIS_WINDOW_BP}bp.tsv"

echo "Wrote: $OUTDIR/cis_raw_lncRNA_PCG_${CIS_WINDOW_BP}bp.tsv"
