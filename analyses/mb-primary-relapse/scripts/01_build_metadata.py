import os
import pandas as pd

# Point these at your local copies of the R2 datagrabber metadata file and the
# raw RNA-seq counts file (not committed to this repo - see ../README.md).
META_FILE = os.environ.get("MB_META_FILE", "ps_avgpres_mbffpeb86_mbffpe_datagrabber.txt")
COUNTS_FILE = os.environ.get("MB_COUNTS_FILE", "MB_primary_relapse.counts.txt")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

fields = ["age", "death", "gender", "group", "histo", "id", "mstage",
          "new.cnv", "new.mut", "pfs", "relapse", "rna_id", "status",
          "subgroup", "tumorid"]

rows = {}
with open(META_FILE) as f:
    for line in f:
        if not line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        tag = parts[0].lstrip("#")
        if tag in fields:
            rows[tag] = parts[2:]  # skip "#tag" and repeated tag col

meta = pd.DataFrame(rows)
meta["sample_id"] = meta["rna_id"].str.upper()
meta = meta.set_index("sample_id")

meta.to_csv(f"{OUT_DIR}/sample_metadata.csv")
print(meta.shape)
print(meta.head(10))
print(meta["group"].value_counts())
print(meta["subgroup"].value_counts())
print(meta["status"].value_counts())

# sanity check against counts file header
with open(COUNTS_FILE) as f:
    header = f.readline().rstrip("\n").split("\t")
counts_ids = set(header)
meta_ids = set(meta.index)
print("in counts not meta:", counts_ids - meta_ids)
print("in meta not counts:", meta_ids - counts_ids)
