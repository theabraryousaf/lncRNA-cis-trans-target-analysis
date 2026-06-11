#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Annotate input genes with GENCODE and create lncRNA/PCG BED files.")
    p.add_argument("--input", required=True, help="Input CSV with gene_id and normalized counts.")
    p.add_argument("--gtf", required=True, help="GENCODE GTF file.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--sample-columns", required=True, help="Comma-separated sample columns.")
    p.add_argument("--lnc-biotypes", default="lncRNA", help="Comma-separated lncRNA gene biotypes. Default: lncRNA")
    p.add_argument("--pcg-biotypes", default="protein_coding", help="Comma-separated PCG biotypes. Default: protein_coding")
    return p.parse_args()


def parse_attrs(attr):
    return dict(re.findall(r'(\S+) "([^"]*)"', attr))


def write_bed(df, path, source_label):
    cols = ["chrom", "start0", "end", "gene_id_base", "gene_symbol", "strand", "gene_type"]
    bed = df[cols].copy()
    bed["source"] = source_label
    bed = bed.dropna(subset=["chrom", "start0", "end"])
    bed["start0"] = bed["start0"].astype(int)
    bed["end"] = bed["end"].astype(int)
    bed = bed.sort_values(["chrom", "start0", "end", "gene_id_base"])
    bed.to_csv(path, sep="\t", header=False, index=False)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample_cols = [x.strip() for x in args.sample_columns.split(",") if x.strip()]
    lnc_biotypes = {x.strip() for x in args.lnc_biotypes.split(",") if x.strip()}
    pcg_biotypes = {x.strip() for x in args.pcg_biotypes.split(",") if x.strip()}

    counts = pd.read_csv(args.input)
    counts["gene_id_base"] = counts["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

    for c in sample_cols:
        counts[c] = pd.to_numeric(counts[c], errors="coerce")

    counts = counts.groupby("gene_id_base", as_index=False)[sample_cols].mean()
    input_ids = set(counts["gene_id_base"])

    gene_records = []

    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs_raw = fields
            if feature != "gene":
                continue

            attrs = parse_attrs(attrs_raw)
            gene_id_versioned = attrs.get("gene_id", "")
            gene_id_base = gene_id_versioned.split(".")[0]

            if gene_id_base not in input_ids:
                continue

            gene_name = attrs.get("gene_name", "")
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))

            gene_records.append({
                "gene_id_base": gene_id_base,
                "gene_id_versioned_gtf": gene_id_versioned,
                "gene_symbol": gene_name if gene_name else gene_id_base,
                "gene_type": gene_type,
                "chrom": chrom,
                "start0": int(start) - 1,
                "end": int(end),
                "strand": strand
            })

    ann = pd.DataFrame(gene_records)
    annotated = counts.merge(ann, on="gene_id_base", how="left")
    annotated["gene_symbol"] = annotated["gene_symbol"].fillna(annotated["gene_id_base"])

    annotated.to_csv(outdir / "significant_genes_annotated.tsv", sep="\t", index=False)
    annotated.to_csv(outdir / "significant_genes_annotated.csv", index=False)

    missing = annotated[annotated["gene_type"].isna()].copy()
    missing.to_csv(outdir / "genes_missing_from_GENCODE.tsv", sep="\t", index=False)

    lnc = annotated[annotated["gene_type"].isin(lnc_biotypes)].copy()
    pcg = annotated[annotated["gene_type"].isin(pcg_biotypes)].copy()

    write_bed(lnc, outdir / "significant_lncRNAs.bed", "input_lncRNA")
    write_bed(pcg, outdir / "significant_protein_coding.bed", "input_protein_coding")

    lnc[["gene_id_base", "gene_symbol", "gene_type"]].to_csv(outdir / "significant_lncRNAs.tsv", sep="\t", index=False)
    pcg[["gene_id_base", "gene_symbol", "gene_type"]].to_csv(outdir / "significant_protein_coding.tsv", sep="\t", index=False)

    print("Input genes:", counts.shape[0])
    print("Annotated genes:", annotated["gene_type"].notna().sum())
    print("Missing from GENCODE:", missing.shape[0])
    print("lncRNAs:", lnc.shape[0])
    print("protein-coding genes:", pcg.shape[0])


if __name__ == "__main__":
    main()
