#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

results = Path("results")


def read_table(path, cols=None):
    path = Path(path)
    if path.exists() and path.stat().st_size > 0:
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame(columns=cols or ["lncRNA_gene_id", "lncRNA_symbol"])


cis = read_table(results / "cis_lncRNA_priority.tsv")
trans = read_table(results / "trans_lncRNA_priority.tsv")

cis = cis.rename(columns={
    "cis_target_count": "cis_target_count_100kb",
    "cis_targets": "cis_targets_100kb"
})
trans = trans.rename(columns={
    "trans_target_count": "trans_target_count_selected",
    "max_abs_PCC": "max_abs_PCC_selected",
    "trans_targets": "trans_targets_selected"
})

final = pd.merge(cis, trans, on=["lncRNA_gene_id", "lncRNA_symbol"], how="outer")

for col in [
    "cis_target_count_100kb",
    "nearest_cis_distance_bp",
    "trans_target_count_selected",
    "positive_trans_target_count",
    "negative_trans_target_count",
    "max_abs_PCC_selected"
]:
    if col not in final.columns:
        final[col] = 0
    final[col] = pd.to_numeric(final[col], errors="coerce").fillna(0)

final["has_cis"] = final["cis_target_count_100kb"] > 0
final["has_trans"] = final["trans_target_count_selected"] > 0


def evidence_type(row):
    if row["has_cis"] and row["has_trans"]:
        return "cis_and_trans"
    if row["has_trans"]:
        return "trans_only"
    if row["has_cis"]:
        return "cis_only"
    return "none"


final["evidence_type"] = final.apply(evidence_type, axis=1)
rank = {"cis_and_trans": 1, "trans_only": 2, "cis_only": 3, "none": 4}
final["evidence_rank"] = final["evidence_type"].map(rank).fillna(99)

final = final.sort_values(
    ["evidence_rank", "trans_target_count_selected", "cis_target_count_100kb", "max_abs_PCC_selected", "nearest_cis_distance_bp"],
    ascending=[True, False, False, False, True]
).drop(columns=["evidence_rank"])

final.to_csv(results / "FINAL_lncRNA_cis_trans_priority.tsv", sep="\t", index=False)
final.to_csv(results / "FINAL_lncRNA_cis_trans_priority.csv", index=False)

print("Final prioritized lncRNAs:", final.shape[0])
print("Evidence type counts:")
print(final["evidence_type"].value_counts().to_string())
print("\nTop candidates:")
show_cols = [
    "lncRNA_gene_id", "lncRNA_symbol", "evidence_type",
    "cis_target_count_100kb", "trans_target_count_selected",
    "positive_trans_target_count", "negative_trans_target_count",
    "max_abs_PCC_selected"
]
print(final[show_cols].head(30).to_string(index=False))
