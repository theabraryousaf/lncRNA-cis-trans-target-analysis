#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

results = Path("results")


def read_table(path):
    path = Path(path)
    if path.exists() and path.stat().st_size > 0:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


ann = read_table(results / "significant_genes_annotated.tsv")
cis_files = sorted(results.glob("cis_targets_*bp_all.tsv"))
cis = read_table(cis_files[0]) if cis_files else pd.DataFrame()
trans_all = read_table(results / "trans_targets_Pearson_all.tsv")
trans_selected = read_table(results / "trans_targets_Pearson_selected.tsv")
final = read_table(results / "FINAL_lncRNA_cis_trans_priority.tsv")

summary = results / "FINAL_cis_trans_analysis_summary.txt"

with open(summary, "w") as out:
    out.write("lncRNA cis/trans target analysis summary\n")
    out.write("=====================================\n\n")

    if not ann.empty:
        out.write(f"Total input genes annotated: {ann.shape[0]}\n")
        out.write(f"lncRNAs: {(ann['gene_type'] == 'lncRNA').sum()}\n")
        out.write(f"protein-coding genes: {(ann['gene_type'] == 'protein_coding').sum()}\n\n")

    out.write(f"Cis pairs within 100 kb: {cis.shape[0]}\n")
    out.write(f"Unique cis target genes: {pd.read_csv(results / 'all_cis_target_genes.tsv', sep='\t').shape[0] if (results / 'all_cis_target_genes.tsv').exists() else 0}\n\n")

    out.write(f"All trans pairs tested: {trans_all.shape[0]}\n")
    out.write(f"Selected trans pairs: {trans_selected.shape[0]}\n")
    if not trans_selected.empty and 'correlation_direction' in trans_selected.columns:
        out.write(f"Positive selected trans pairs: {(trans_selected['correlation_direction'] == 'positive').sum()}\n")
        out.write(f"Negative selected trans pairs: {(trans_selected['correlation_direction'] == 'negative').sum()}\n")
    out.write(f"Unique trans target genes: {pd.read_csv(results / 'all_trans_target_genes.tsv', sep='\t').shape[0] if (results / 'all_trans_target_genes.tsv').exists() else 0}\n\n")

    if not final.empty:
        out.write("Evidence type counts:\n")
        out.write(final["evidence_type"].value_counts().to_string())
        out.write("\n\nTop 30 prioritized lncRNAs:\n")
        cols = [
            "lncRNA_gene_id", "lncRNA_symbol", "evidence_type",
            "cis_target_count_100kb", "trans_target_count_selected",
            "positive_trans_target_count", "negative_trans_target_count",
            "max_abs_PCC_selected"
        ]
        cols = [c for c in cols if c in final.columns]
        out.write(final[cols].head(30).to_string(index=False))
        out.write("\n")

print("Wrote:", summary)
print(summary.read_text())
