#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Create log2 normalized expression matrix for Pearson correlation.")
    p.add_argument("--input", required=True, help="Input CSV.")
    p.add_argument("--sample-columns", required=True, help="Comma-separated sample columns.")
    p.add_argument("--output", required=True, help="Output TSV.")
    return p.parse_args()


def main():
    args = parse_args()
    sample_cols = [x.strip() for x in args.sample_columns.split(",") if x.strip()]
    df = pd.read_csv(args.input)

    if "gene_id" not in df.columns:
        raise SystemExit("ERROR: Input file must contain gene_id column.")

    df["gene_id"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if (df[sample_cols] < 0).any().any():
        raise SystemExit("ERROR: Negative expression values found. log2(count + 1) requires non-negative values.")

    df = df.groupby("gene_id", as_index=False)[sample_cols].mean()
    df[sample_cols] = np.log2(df[sample_cols] + 1)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    print("Wrote:", out)
    print("Genes:", df.shape[0])
    print("Samples:", len(sample_cols))


if __name__ == "__main__":
    main()
