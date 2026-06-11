#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Validate input normalized-count CSV for lncRNA cis/trans analysis.")
    p.add_argument("--input", required=True, help="Input CSV with gene_id and sample columns.")
    p.add_argument("--sample-columns", required=True, help="Comma-separated sample columns.")
    return p.parse_args()


def main():
    args = parse_args()
    infile = Path(args.input)
    sample_cols = [x.strip() for x in args.sample_columns.split(",") if x.strip()]

    if not infile.exists():
        raise SystemExit(f"ERROR: Missing input file: {infile}")

    df = pd.read_csv(infile)

    if "gene_id" not in df.columns:
        raise SystemExit("ERROR: File must contain a column named gene_id.")

    if len(sample_cols) < 3:
        raise SystemExit("ERROR: At least 3 sample columns are required for Pearson correlation.")

    missing = [c for c in sample_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: Missing sample columns: {missing}")

    df["gene_id_base"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    duplicated = df["gene_id_base"].duplicated().sum()
    if duplicated:
        print(f"WARNING: {duplicated} duplicated Ensembl base IDs found. They will be averaged later.")

    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if df[sample_cols].isna().any().any():
        bad_cols = df[sample_cols].columns[df[sample_cols].isna().any()].tolist()
        raise SystemExit(f"ERROR: Non-numeric or missing values found in sample columns: {bad_cols}")

    if (df[sample_cols] < 0).any().any():
        raise SystemExit("ERROR: Negative expression values found. log2(count + 1) requires non-negative values.")

    print("Input file OK")
    print("Genes:", df.shape[0])
    print("Samples:", len(sample_cols))
    print("Sample columns:", ", ".join(sample_cols))


if __name__ == "__main__":
    main()
