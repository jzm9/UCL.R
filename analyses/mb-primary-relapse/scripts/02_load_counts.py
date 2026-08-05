import os
import pandas as pd
import numpy as np

# Point this at your local copy of the raw RNA-seq counts file (not committed
# to this repo - see ../README.md).
COUNTS_FILE = os.environ.get("MB_COUNTS_FILE", "MB_primary_relapse.counts.txt")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

counts = pd.read_csv(COUNTS_FILE, sep="\t", index_col=0)
print("counts shape:", counts.shape)

genes_of_interest = ["LRG1", "CD74", "MIF", "EMILIN3", "ENG", "ITGB1", "PTK2"]
for g in genes_of_interest:
    print(g, "present:", g in counts.index)

# CPM + log2 normalization (library-size normalized, standard for raw RNA-seq counts)
libsize = counts.sum(axis=0)
cpm = counts.div(libsize, axis=1) * 1e6
logcpm = np.log2(cpm + 1)

logcpm.to_csv(f"{OUT_DIR}/logcpm_matrix.csv")
print("logcpm shape:", logcpm.shape)
print(logcpm.loc[genes_of_interest].T.describe())
