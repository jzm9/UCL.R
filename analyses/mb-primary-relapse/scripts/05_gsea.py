import os
import pandas as pd
import numpy as np
from scipy import stats
import gseapy as gp

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

meta = pd.read_csv(f"{OUT_DIR}/sample_metadata.csv", index_col=0)
logcpm = pd.read_csv(f"{OUT_DIR}/logcpm_matrix.csv", index_col=0)

# gene sets: verified MSigDB Hallmark v2023.1 (SHA-confirmed across 4 independent GitHub
# mirrors of the official Broad Institute release) + a lower-confidence supplementary
# curated panel + the user's specific 7-gene panel.
HALLMARK_ANGIOGENESIS = ["VCAN", "POSTN", "FSTL1", "LRPAP1", "STC1", "LPL", "VEGFA", "PF4",
    "THBD", "FGFR1", "TNFRSF21", "CCND2", "COL5A2", "ITGAV", "SERPINA5", "KCNJ8", "APP",
    "JAG1", "COL3A1", "SPP1", "NRP1", "OLR1", "PDGFA", "PTK2", "SLCO2A1", "PGLYRP1", "VAV2",
    "S100A4", "MSX1", "VTN", "TIMP1", "APOH", "PRG2", "JAG2", "LUM", "CXCL6"]

CURATED_ANGIOGENESIS_SUPPLEMENTARY = ["VEGFA", "VEGFB", "VEGFC", "VEGFD", "KDR", "FLT1",
    "FLT4", "NRP1", "NRP2", "PDGFA", "PDGFB", "FGF1", "FGF2", "ANGPT1", "ANGPT2", "TEK",
    "TIE1", "THBS1", "SERPINE1", "MMP2", "MMP9", "TIMP1", "TIMP2", "TIMP3", "COL18A1",
    "COL4A2", "SPP1", "VTN", "ITGAV", "ITGB3", "CDH5", "PECAM1", "MCAM", "VCAM1", "SELP",
    "CXCL8", "CXCL12", "CCL2", "LPL", "PLAU", "PLAUR", "NOS3", "HIF1A", "EPAS1", "ACVRL1",
    "ENG", "TGFBR1", "SMAD4", "EFNB2", "EPHB4", "ROBO4"]

USER_PANEL = ["LRG1", "CD74", "MIF", "EMILIN3", "ENG", "ITGB1", "PTK2"]

gene_sets = {
    "HALLMARK_ANGIOGENESIS": HALLMARK_ANGIOGENESIS,
    "CURATED_ANGIOGENESIS_SUPPLEMENTARY": CURATED_ANGIOGENESIS_SUPPLEMENTARY,
    "USER_ANGIOGENESIS_PANEL": USER_PANEL,
}

group_order = ["group_3", "group_4", "shh"]
group_label = {"group_3": "Group 3 MB", "group_4": "Group 4 MB", "shh": "SHH-MB"}

# filter to expressed genes (mean CPM > 1 i.e. mean log2CPM > log2(2)) for stable ranking
expressed = logcpm.index[logcpm.mean(axis=1) > 1]
logcpm_f = logcpm.loc[expressed]
print(f"genes after expression filter: {logcpm_f.shape[0]} / {logcpm.shape[0]}")

# Pair by tumorid globally, then bucket each pair by the PRIMARY sample's group.
# Two Group 4->Group 3 switch-at-relapse pairs exist (matches the paper's report of
# ~5% of Group 3/4 cases switching molecular group at relapse) - bucketing by primary
# group avoids splitting/dropping those pairs.
prim_all = meta[meta["status"] == "prim"]
recur_all = meta[meta["status"] == "recur"]
pair_group = prim_all.set_index("tumorid")["group"]

all_results = []
for grp in group_order:
    tumorids = pair_group[pair_group == grp].index
    prim_sub = prim_all[prim_all["tumorid"].isin(tumorids)].sort_values("tumorid")
    recur_sub = recur_all[recur_all["tumorid"].isin(tumorids)].sort_values("tumorid")
    prim_ids = prim_sub.index
    recur_ids = recur_sub.index
    assert (prim_sub["tumorid"].values == recur_sub["tumorid"].values).all()

    p = logcpm_f[prim_ids].values
    r = logcpm_f[recur_ids].values
    diff = r - p
    n = diff.shape[1]
    mean_diff = diff.mean(axis=1)
    sd_diff = diff.std(axis=1, ddof=1)
    sd_diff[sd_diff == 0] = np.nan
    tstat = mean_diff / (sd_diff / np.sqrt(n))
    rnk = pd.Series(tstat, index=logcpm_f.index).dropna().sort_values(ascending=False)
    rnk_path = f"{OUT_DIR}/rank_{grp}.rnk"
    rnk.to_csv(rnk_path, sep="\t", header=False)

    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=3,
        max_size=2000,
        permutation_num=1000,
        outdir=None,
        seed=42,
        verbose=False,
    )
    res = pre_res.res2d.copy()
    res["group"] = group_label[grp]
    all_results.append(res)
    print(f"\n=== {group_label[grp]} (n={n} pairs, ranked by paired t-stat, relapse vs primary) ===")
    print(res[["Term", "NES", "NOM p-val", "FDR q-val", "Tag %"]].to_string(index=False))

final = pd.concat(all_results, ignore_index=True)
final.to_csv(f"{OUT_DIR}/gsea_angiogenesis_results.csv", index=False)
print("\nSaved:", f"{OUT_DIR}/gsea_angiogenesis_results.csv")
