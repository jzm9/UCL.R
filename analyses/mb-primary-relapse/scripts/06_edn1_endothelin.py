import os
import numpy as np
import pandas as pd
from scipy import stats

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

meta = pd.read_csv(f"{OUT_DIR}/sample_metadata.csv", index_col=0)
logcpm = pd.read_csv(f"{OUT_DIR}/logcpm_matrix.csv", index_col=0)

# Endothelial marker score (PECAM1 is absent from this count matrix's annotation)
ENDOTHELIAL = ["CDH5", "KDR", "VWF", "CLDN5", "ESAM"]
ez = logcpm.loc[ENDOTHELIAL]
ez = ez.sub(ez.mean(axis=1), axis=0).div(ez.std(axis=1), axis=0)
endo_score = ez.mean()

features = {g: logcpm.loc[g] for g in ["EDN1", "EDNRA", "EDNRB", "HIF1A"]}
features["endothelial_score"] = endo_score
# EDN1 relative to vessel content: residual after regressing EDN1 on the endothelial score
fit = stats.linregress(endo_score, logcpm.loc["EDN1"])
features["EDN1_adj_endothelial"] = logcpm.loc["EDN1"] - (fit.intercept + fit.slope * endo_score)

# Pair by tumorid, bucket each pair by its primary's group (same as 05_gsea.py)
prim = meta[meta["status"] == "prim"]
recur = meta[meta["status"] == "recur"]
pair_group = prim.set_index("tumorid")["group"]
group_label = {"group_3": "Group 3 MB", "group_4": "Group 4 MB", "shh": "SHH-MB"}

records = []
for grp in ["group_3", "group_4", "shh"]:
    tumorids = pair_group[pair_group == grp].index
    p = prim[prim["tumorid"].isin(tumorids)].sort_values("tumorid")
    r = recur[recur["tumorid"].isin(tumorids)].sort_values("tumorid")
    assert (p["tumorid"].values == r["tumorid"].values).all()
    for name, values in features.items():
        diff = values[r.index].values - values[p.index].values
        records.append({
            "group": group_label[grp],
            "feature": name,
            "n_pairs": len(diff),
            "median_delta_relapse_minus_primary": np.median(diff),
            "n_up": int((diff > 0).sum()),
            "n_down": int((diff < 0).sum()),
            "wilcoxon_p": stats.wilcoxon(diff).pvalue,
        })

paired = pd.DataFrame(records)
paired.to_csv(f"{OUT_DIR}/edn1_paired_by_subgroup.csv", index=False)
print(paired.to_string(index=False))

# Cross-sectional correlation of EDN1 with vascular/hypoxia signals, all 86 samples
corr = []
for name, v in [("endothelial_score", endo_score)] + [(g, logcpm.loc[g]) for g in
                ["CDH5", "VWF", "VEGFA", "HIF1A", "EDNRA", "EDNRB", "LRG1"]]:
    rho, pval = stats.spearmanr(logcpm.loc["EDN1"], v)
    corr.append({"feature": name, "spearman_rho": rho, "p": pval})
corr = pd.DataFrame(corr)
corr.to_csv(f"{OUT_DIR}/edn1_correlations.csv", index=False)
print(corr.to_string(index=False))
