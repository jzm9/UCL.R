import os
import pandas as pd
import numpy as np
from scipy import stats
import json

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

meta = pd.read_csv(f"{OUT_DIR}/sample_metadata.csv", index_col=0)
logcpm = pd.read_csv(f"{OUT_DIR}/logcpm_matrix.csv", index_col=0)

genes = ["LRG1", "CD74", "MIF", "EMILIN3", "ENG", "ITGB1", "PTK2"]
gene_label = {"PTK2": "PTK2 (FAK)", "ENG": "ENG (Endoglin)"}

group_map = {"shh": "SHH-MB", "group_3": "Group 3 MB", "group_4": "Group 4 MB"}
group_order = ["group_3", "group_4", "shh"]  # Group 3 first (most important)

expr = logcpm.loc[genes].T  # samples x genes
expr = expr.join(meta[["group", "status", "tumorid"]])

# z-score each gene across the whole cohort (for the composite angiogenesis score)
z = (logcpm.loc[genes].T - logcpm.loc[genes].T.mean()) / logcpm.loc[genes].T.std()
expr["angio_score"] = z.mean(axis=1)

records = []
paired_data = {}  # for plotting: group -> gene -> list of (tumorid, prim_val, recur_val)

for grp in group_order:
    sub = expr[expr["group"] == grp]
    pivot_prim = sub[sub["status"] == "prim"].set_index("tumorid")
    pivot_recur = sub[sub["status"] == "recur"].set_index("tumorid")
    common = pivot_prim.index.intersection(pivot_recur.index)
    n_pairs = len(common)

    for gene in genes + ["angio_score"]:
        p_vals = pivot_prim.loc[common, gene].values
        r_vals = pivot_recur.loc[common, gene].values
        diff = r_vals - p_vals
        if n_pairs >= 2 and np.any(diff != 0):
            wstat, wp = stats.wilcoxon(p_vals, r_vals)
        else:
            wp = np.nan
        tstat, tp = stats.ttest_rel(p_vals, r_vals) if n_pairs >= 2 else (np.nan, np.nan)

        records.append({
            "group": group_map[grp],
            "gene": gene,
            "n_pairs": n_pairs,
            "median_primary": np.median(p_vals),
            "median_relapse": np.median(r_vals),
            "log2FC_relapse_vs_primary": np.median(diff),
            "n_increased_at_relapse": int(np.sum(diff > 0)),
            "n_decreased_at_relapse": int(np.sum(diff < 0)),
            "paired_ttest_p": tp,
            "wilcoxon_p": wp,
        })

        paired_data.setdefault(grp, {})[gene] = {
            "tumorid": list(common),
            "primary": p_vals.tolist(),
            "relapse": r_vals.tolist(),
        }

results = pd.DataFrame(records)
results.to_csv(f"{OUT_DIR}/gene_stats_by_subgroup.csv", index=False)
pd.set_option("display.width", 140)
print(results.to_string(index=False))

with open(f"{OUT_DIR}/paired_data.json", "w") as f:
    json.dump(paired_data, f)
