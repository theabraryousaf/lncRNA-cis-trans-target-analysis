#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Parse BEDTools cis window output.")
    p.add_argument("--raw", required=True, help="Raw bedtools window output.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--window-bp", required=True, type=int, help="Cis window size used.")
    return p.parse_args()


def clean_symbol(x, gene_id):
    if pd.isna(x) or str(x).strip().lower() in {"", "nan", "none"}:
        return gene_id
    return str(x).strip()


def main():
    args = parse_args()
    raw_file = Path(args.raw)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cols = [
        "lnc_chrom", "lnc_start0", "lnc_end", "lncRNA_gene_id", "lncRNA_symbol",
        "lnc_strand", "lnc_gene_type", "lnc_source",
        "target_chrom", "target_start0", "target_end", "target_gene_id", "target_symbol",
        "target_strand", "target_gene_type", "target_source"
    ]

    pair_file = outdir / f"cis_targets_{args.window_bp}bp_all.tsv"
    priority_file = outdir / "cis_lncRNA_priority.tsv"
    gene_file = outdir / "all_cis_target_genes.tsv"
    gene_txt = outdir / "all_cis_target_genes.txt"

    if not raw_file.exists() or raw_file.stat().st_size == 0:
        empty_pairs = pd.DataFrame(columns=cols + ["distance_bp", "target_position_relative_to_lncRNA"])
        empty_pairs.to_csv(pair_file, sep="\t", index=False)
        empty_pairs.to_csv(pair_file.with_suffix(".csv"), index=False)
        pd.DataFrame(columns=["lncRNA_gene_id", "lncRNA_symbol", "cis_target_count", "nearest_cis_distance_bp", "cis_targets"]).to_csv(priority_file, sep="\t", index=False)
        pd.DataFrame(columns=["target_gene_id", "target_symbol"]).to_csv(gene_file, sep="\t", index=False)
        gene_txt.write_text("")
        print("No cis pairs found.")
        return

    raw = pd.read_csv(raw_file, sep="\t", header=None, names=cols)

    raw["lncRNA_symbol"] = [clean_symbol(s, g) for s, g in zip(raw["lncRNA_symbol"], raw["lncRNA_gene_id"])]
    raw["target_symbol"] = [clean_symbol(s, g) for s, g in zip(raw["target_symbol"], raw["target_gene_id"])]

    def distance(row):
        if row["target_end"] <= row["lnc_start0"]:
            return int(row["lnc_start0"] - row["target_end"])
        if row["lnc_end"] <= row["target_start0"]:
            return int(row["target_start0"] - row["lnc_end"])
        return 0

    def relative_position(row):
        if row["target_end"] <= row["lnc_start0"]:
            side = "left"
        elif row["lnc_end"] <= row["target_start0"]:
            side = "right"
        else:
            return "overlapping"

        if row["lnc_strand"] == "+":
            return "upstream" if side == "left" else "downstream"
        if row["lnc_strand"] == "-":
            return "downstream" if side == "left" else "upstream"
        return side

    raw["distance_bp"] = raw.apply(distance, axis=1)
    raw["target_position_relative_to_lncRNA"] = raw.apply(relative_position, axis=1)

    raw = raw.drop_duplicates(["lncRNA_gene_id", "target_gene_id", "distance_bp", "target_position_relative_to_lncRNA"])
    raw = raw.sort_values(["lncRNA_symbol", "distance_bp", "target_symbol"])

    raw.to_csv(pair_file, sep="\t", index=False)
    raw.to_csv(pair_file.with_suffix(".csv"), index=False)

    priority = (
        raw.groupby(["lncRNA_gene_id", "lncRNA_symbol"], dropna=False)
        .agg(
            cis_target_count=("target_gene_id", "nunique"),
            nearest_cis_distance_bp=("distance_bp", "min"),
            cis_targets=("target_symbol", lambda x: ";".join(sorted(set(map(str, x)))))
        )
        .reset_index()
        .sort_values(["cis_target_count", "nearest_cis_distance_bp"], ascending=[False, True])
    )

    priority.to_csv(priority_file, sep="\t", index=False)
    priority.to_csv(priority_file.with_suffix(".csv"), index=False)

    genes = raw[["target_gene_id", "target_symbol"]].drop_duplicates().sort_values("target_symbol")
    genes.to_csv(gene_file, sep="\t", index=False)
    genes.to_csv(gene_file.with_suffix(".csv"), index=False)
    gene_txt.write_text("\n".join(genes["target_symbol"].astype(str)) + "\n")

    print("Cis lncRNA-PCG pairs:", raw.shape[0])
    print("lncRNAs with cis targets:", priority.shape[0])
    print("Unique cis target genes:", genes.shape[0])


if __name__ == "__main__":
    main()
