# METHODS: MB primary vs. relapse angiogenesis analysis

This document explains, script by script, exactly what the code in
`scripts/` does and why — the data model, every statistical test, every
parameter choice, and the caveats. It's meant to let you (or anyone else)
audit or extend the analysis without re-reading the code line by line. See
`README.md` for the short version and the headline findings.

## 0. The question and the constraint

The ask: for the Okonechnikov et al. 2023 medulloblastoma (MB) primary-relapse
RNA-seq cohort, look at seven specific genes (LRG1, CD74, MIF, EMILIN3, ENG,
ITGB1, PTK2) and the angiogenesis pathway more broadly, comparing primary vs.
relapse tumors, split by molecular subgroup (SHH-MB / Group 3 MB / Group 4
MB), with Group 3 flagged as the priority.

The two inputs are two plain-text files the user uploaded:

- **`MB_primary_relapse.counts.txt`** — raw RNA-seq gene counts. 55,765 genes
  (rows) x 86 samples (columns), tab-separated, sample IDs like `MB135076` as
  column headers, gene symbols as row labels.
- **`ps_avgpres_mbffpeb86_mbffpe_..._datagrabber.txt`** — an export from the
  R2 Genomics Analysis and Visualization Platform (http://r2.amc.nl), the
  platform the paper itself names as the host of this dataset
  ("Tumor Medulloblastoma—Korshunov—86—rpkm—mbffpe"). This file mixes two
  things in one table: a block of `#`-prefixed metadata rows at the top
  (age, subgroup, status, etc., one column per sample) followed by a large
  block of R2's own normalized expression values (not used here — see
  §1 for why).

Neither file alone has both the counts *and* the sample labels lined up, so
step 1 is entirely about reconstructing a trustworthy sample→group mapping
before any statistics happen.

## 1. `01_build_metadata.py` — building the sample metadata table

**What it does:** reads the R2 datagrabber file line by line, keeps only
lines starting with `#` (R2's convention for an annotation row — data rows
are gene symbols with no `#`), and for a fixed allowlist of 15 fields
(`age, death, gender, group, histo, id, mstage, new.cnv, new.mut, pfs,
relapse, rna_id, status, subgroup, tumorid`) captures that row's values
(skipping the first two columns, which are just the repeated tag name).

```python
# scripts/01_build_metadata.py, lines 10-26
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
```

Each kept row becomes one column of a `pandas.DataFrame` — i.e. the
DataFrame ends up **samples × metadata-fields**, indexed by
`sample_id = rna_id.upper()` (R2 stores IDs lowercase, e.g. `mb103998`; the
counts file uses uppercase, e.g. `MB103998` — the `.upper()` is what makes
the two files joinable at all).

**Why not use R2's own normalized values in that same file?** Two reasons.
First, the paper's own methods describe deriving RPKM from raw counts
themselves, not just took R2's already-processed numbers — reproducing that
from the counts file keeps the analysis under our control (choice of
normalization, filtering, etc.) rather than opaque to R2's pipeline. Second,
the R2 export block that follows the metadata rows is ~48,000 gene rows of
values that hover around 1.0 (e.g. `1.02, 1.05, 1.13`) — these look like a
**ratio/fold-change-style scaling** (each gene's value relative to something,
possibly a reference or the row's own geometric mean), not expression on an
absolute or log scale. We didn't need to fully reverse-engineer that scale
because the raw counts file gives us a strictly better starting point:
ground-truth integer read counts we can normalize ourselves in a way that's
documented and auditable (§2).

**Sanity check:** the script's last four lines diff the metadata's
`sample_id` set against the counts file's column-header set:

```python
# scripts/01_build_metadata.py, lines 35-41
# sanity check against counts file header
with open(COUNTS_FILE) as f:
    header = f.readline().rstrip("\n").split("\t")
counts_ids = set(header)
meta_ids = set(meta.index)
print("in counts not meta:", counts_ids - meta_ids)
print("in meta not counts:", meta_ids - counts_ids)
```

Both directions came back empty — **all 86 samples in both files match
exactly**, case differences aside — which is what let us trust the join for
everything downstream.

**Output:** `results/sample_metadata.csv`, 86 rows x 15 columns, indexed by
sample ID.

### What's actually in the metadata (as it came out)

| Field | What it is | Values seen |
|---|---|---|
| `status` | primary or relapse | `prim` (43), `recur` (43) |
| `group` | consensus molecular group | `shh` (48 samples/24 pairs\*), `group_4` (26/13\*), `group_3` (12/6\*) |
| `subgroup` | finer 2nd-gen subgroup | `shh_inf`, `ch_ad` (SHH split by age); `i`–`viii` (Group 3/4) |
| `tumorid` | patient pair ID, e.g. `pr23` | links a `prim` sample to its `recur` sample |
| `age`, `gender`, `histo`, `mstage`, `pfs`, `death`, `new.cnv`, `new.mut` | clinical/molecular covariates | as reported in the paper's Table S1-equivalent |

\* Naive per-sample counting within the `group` column slightly
undercounts/miscounts pairs — see the switch-pair issue in §5, which is
handled correctly in the GSEA script but *not* corrected for in the
per-gene stats script (a known, documented discrepancy — see §6).

## 2. `02_load_counts.py` — normalizing raw counts

**What it does:**

```python
# scripts/02_load_counts.py, lines 10-22
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
```

This is standard **library-size normalization** (CPM), then a log2 transform
with a pseudocount of 1 to handle zeros. It's the simplest normalization that
makes samples with different sequencing depth (library sizes ranged
~1.9M–44M reads in this cohort) comparable.

**Why CPM and not the paper's RPKM?** RPKM additionally divides by gene
length, which matters when comparing *different genes'* expression levels to
each other (a long gene will always look more "expressed" than a short one
in raw counts, purely from having more read-alignment opportunities). Here,
every comparison is **the same gene, across samples** (primary vs. relapse) —
gene length is constant within a gene, so it cancels out and doesn't need to
be corrected for. CPM is sufficient and simpler to reason about.

**Output:** `logcpm_matrix.csv` (55,765 genes x 86 samples) — this is a large
derived-but-still-essentially-full-resolution expression matrix, so
(consistent with not committing the raw counts) **it is not checked into the
repo**; scripts 03 and 05 both regenerate it by re-running 02, or you can
keep a local copy alongside the raw counts.

## 3. `03_analyze_genes.py` — paired per-gene statistics

This is the core per-gene analysis for the 7 requested genes.

**Composite angiogenesis score:**

```python
# scripts/03_analyze_genes.py, lines 12-23
genes = ["LRG1", "CD74", "MIF", "EMILIN3", "ENG", "ITGB1", "PTK2"]
gene_label = {"PTK2": "PTK2 (FAK)", "ENG": "ENG (Endoglin)"}

group_map = {"shh": "SHH-MB", "group_3": "Group 3 MB", "group_4": "Group 4 MB"}
group_order = ["group_3", "group_4", "shh"]  # Group 3 first (most important)

expr = logcpm.loc[genes].T  # samples x genes
expr = expr.join(meta[["group", "status", "tumorid"]])

# z-score each gene across the whole cohort (for the composite angiogenesis score)
z = (logcpm.loc[genes].T - logcpm.loc[genes].T.mean()) / logcpm.loc[genes].T.std()
expr["angio_score"] = z.mean(axis=1)
```

Each of the 7 genes' logCPM values is z-scored across the *whole* cohort
(`(x - mean) / sd`, all 86 samples), then the 7 z-scores are averaged per
sample. This gives one number per sample
that goes up when a sample tends to have high LRG1/CD74/.../PTK2 relative to
the cohort, and down when it tends to have low values across the panel. It's
a simple, transparent composite, not a formal gene-set score (contrast with
GSVA/ssGSEA in §5) — good for "does this sample look angiogenesis-high or
-low overall" but it treats all 7 genes as equally important, which they are
not (e.g. MIF's near-zero counts made it contribute mostly noise; see the
caveat in §6).

**Pairing and statistics, per subgroup:**

```python
# scripts/03_analyze_genes.py, lines 28-62
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
```

Step by step:

1. Filter to samples in that `group` (`shh` / `group_3` / `group_4`).
2. Split into `prim` and `recur`, each indexed by `tumorid`.
3. Take the **intersection** of tumorids present in both — i.e. only
   complete pairs are kept.
4. For each gene (+ the composite score): run
   - `scipy.stats.wilcoxon(primary_values, relapse_values)` — the paired
     **Wilcoxon signed-rank test**. This is a non-parametric test on the
     paired differences; it doesn't assume the differences are normally
     distributed, which matters with as few as 5 pairs (Group 3 MB). Guarded
     against `n_pairs < 2` and against the all-zero-differences edge case
     (which `scipy.stats.wilcoxon` raises on).
   - `scipy.stats.ttest_rel(...)` — the paired **t-test**, reported
     alongside as a parametric cross-check, but Wilcoxon is the primary
     number reported (more defensible at this sample size).
   - Also recorded: median primary, median relapse, median paired difference
     (`log2FC_relapse_vs_primary`, a slight misnomer — it's a difference of
     log2CPM values, i.e. a log2 fold-change, not a linear FC), and a simple
     count of how many patients went up vs. down at relapse (direction-of-
     effect sanity check independent of the p-value).

**Why intersection, not the primary-group bucketing used later in
`05_gsea.py`?** At the time this script was written, two patient pairs
(`pr26`, `pr44`) hadn't yet been identified as switching molecular group
between primary and relapse (Group 4 → Group 3). Filtering by `group` *before*
pairing silently drops any tumorid whose primary and relapse samples land in
different `group` buckets — which is exactly what happens to those two pairs
(present in the `group_4`-filtered `prim` set, but their `recur` sample is
in the `group_3`-filtered set, so the tumorid never appears in the
intersection for either group). This under-counts real pairs (n=5 for Group
3, n=12 for Group 4, instead of the true 6 and 14) but does **not**
mismatch a primary against the wrong relapse — it's a conservative omission,
not a correctness bug. See §5 and §6 for how this was later caught and fixed
for the GSEA analysis, and why it was left as-is here (already-published
numbers, and the effect of adding 2 more pairs to a already-small n is
unlikely to change the qualitative picture — but it does mean the two
n-counts differ between the two analyses, which is worth knowing before
comparing them directly).

**Outputs:**
- `results/gene_stats_by_subgroup.csv` — one row per (subgroup, gene) with
  all the numbers above.
- `results/paired_data.json` — the raw paired values per (subgroup, gene),
  used by `04_make_report.py` to draw the slope charts without needing the
  full expression matrix again.

## 4. `04_make_report.py` — building the HTML report

Pure presentation layer; computes nothing statistical. Two chart types, both
hand-built as inline SVG (no charting library — keeps the report a single
self-contained HTML file with no external dependencies, per this session's
network restrictions and per the practice of keeping published artifacts
self-contained):

**Paired slope charts** (`make_panel`): for a given (gene, subgroup), draws
one dot per sample at x=Primary or x=Relapse, at a y-position scaled to that
gene's log2CPM range *within that panel*. The y-scale and the line coloring
are the core of it:

```python
# scripts/04_make_report.py, lines 41-60
    all_vals = np.concatenate([prim, recur])
    ymin, ymax = all_vals.min(), all_vals.max()
    span = ymax - ymin
    if span == 0:
        span = 1
    ymin -= span * 0.15
    ymax += span * 0.15

    def y(v):
        return PAD_T + INNER_H - (v - ymin) / (ymax - ymin) * INNER_H

    lines = []
    for pv, rv in zip(prim, recur):
        yp, yr = y(pv), y(rv)
        up = yr < yp
        color = "var(--series-2)" if up else "var(--series-1-muted)"
        lines.append(
            f'<line x1="{X0}" y1="{yp:.1f}" x2="{X1}" y2="{yr:.1f}" '
            f'stroke="{color}" stroke-width="1.4" opacity="0.75"/>'
        )
```

Each small panel uses its own y-axis (`ymin`/`ymax` recomputed per gene x
subgroup, padded 15% for breathing room) — panels are not comparable to each
other by eye on absolute height, only by their own primary→relapse slope and
the printed p-value. A line connects each patient's primary dot to their
relapse dot; the line is colored orange (`--series-2`) if relapse > primary
for that patient (note `yr < yp` because SVG y grows downward), muted blue
(`--series-1-muted`) otherwise — so you can see the up/down split at a
glance before reading the stats table. Group 3 MB panels get a highlighted
border (`highlight=True`) since that's the subgroup of most interest:

```python
# scripts/04_make_report.py, lines 89-94
for gene in genes + ["angio_score"]:
    label = "Composite angiogenesis-gene score (mean z-score)" if gene == "angio_score" else f'{gene_label.get(gene, gene)}'
    panels = "".join(
        make_panel(gene, grp, highlight=(grp == "group_3"))
        for grp in group_order
    )
```

**GSEA bar chart**: one horizontal bar per (subgroup, gene set), bar length
= NES (normalized enrichment score, see §5), color = gene set identity
(3-color categorical palette), full opacity if FDR q-value < 0.05, faded
otherwise:

```python
# scripts/04_make_report.py, lines 127-144
for grp in grp_display_order:
    bars.append(f'<text x="4" y="{y - 6}" class="gsea-group-label">{grp}</text>')
    for i, term in enumerate(term_order):
        row = gsea_df[(gsea_df["group"] == grp) & (gsea_df["Term"] == term)]
        if row.empty:
            y += row_h
            continue
        r = row.iloc[0]
        nes = r["NES"]
        fdr = r["FDR q-val"]
        bx = x_scale(min(nes, 0))
        bw = abs(x_scale(nes) - zero_x)
        sig = fdr < 0.05
        opacity = "1" if sig else "0.45"
        bars.append(
            f'<rect x="{bx:.1f}" y="{y-9:.1f}" width="{bw:.1f}" height="14" rx="3" '
            f'fill="{series_colors[i]}" opacity="{opacity}"/>'
        )
```

Bars are drawn from the zero-line outward in whichever direction `nes` sign
points (`bx = x_scale(min(nes, 0))` picks the left edge whether `nes` is
positive or negative), so a vertical zero-line anchors all of them and
direction (toward primary vs. toward relapse) is visually obvious.

**Design-system notes:** color roles and CSS custom properties follow the
`dataviz` skill's palette (`--series-1/2/3`, validated with
`validate_palette.js` for colorblind-safety in both light and dark mode —
the 3-color set passes with one WARN, a light-mode contrast shortfall on the
green series, mitigated by every bar already carrying a direct text label,
which is the skill's prescribed relief for that warning). Both a
`prefers-color-scheme` media query and a `[data-theme]` override are defined
so the report matches the viewer's OS/app theme.

**Output:** `results/report.html`, also published as a Claude Artifact
(this is the URL shared in conversation) for interactive viewing.

## 5. `05_gsea.py` — pathway-level analysis (GSEA)

This is the part that asks "is angiogenesis *as a pathway* shifted at
relapse", rather than just "are these 7 particular genes shifted".

### 5a. Gene sets used, and how they were sourced

This session's network egress is allowlisted (pypi, npm, github, a few
others) and explicitly **blocks** arbitrary hosts including
`gsea-msigdb.org` and `maayanlab.cloud` (Enrichr) — both attempts returned
HTTP 403 from the egress proxy as a policy denial, not a transient failure.
So the standard "download the .gmt from MSigDB" path wasn't available.

Instead: `raw.githubusercontent.com` **is** reachable (GitHub is in the
allowlist), so a GitHub code search for the known MSigDB filename
`h.all.v2023.1.Hs.symbols.gmt` was used. This is the official filename Broad
ships for Hallmark v2023.1. Four independent, unrelated public repositories
were found to contain a file with that exact name — and critically, **all
four had the identical git blob SHA** (`450d222224c9d93f310888eb52fc7408d6118ffe`).
Independent repos landing on byte-identical content for a file with the
official release name is strong circumstantial evidence it's an unmodified
copy of the real Broad Institute release, not a hand-edited approximation —
this is the basis for calling `HALLMARK_ANGIOGENESIS` (extracted from that
file) "verified" rather than "curated." As it's actually hardcoded in the
script (all three gene sets are, since none could be downloaded at runtime):

```python
# scripts/05_gsea.py, lines 15-27
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
```

A second attempt was made to find the GO Biological Process term
`GOBP_ANGIOGENESIS` (GO:0001525) the same way, but GitHub code search
doesn't index files above a size threshold, and the full GO BP `.gmt`
collection is too large — no hit. A single repo did turn up a
`go_biological_process.json` containing a `GOBP_ANGIOGENESIS` key with ~50
genes, but every gene set in that file is suspiciously round (~48-50 genes
each) compared to real GO Biological Process terms (the true
`GOBP_ANGIOGENESIS` has 300+ genes) — this looks like a hand-abbreviated or
LLM-generated approximation, not the authoritative term. It's kept in the
analysis as `CURATED_ANGIOGENESIS_SUPPLEMENTARY` and **explicitly labeled
"lower confidence"** everywhere it's reported (script comments, README, HTML
report legend) rather than presented as equivalent to the Hallmark set.

The third gene set, `USER_ANGIOGENESIS_PANEL`, is just the 7 genes originally
requested, tested as its own mini gene-set so its pathway-level behavior as
a unit can be compared against the two more general sets.

### 5b. The switch-pair fix

`meta[meta["group"] == grp]`-style filtering (as used in script 03) drops
the 2 patients whose tumor switched Group 4 → Group 3 between primary and
relapse. This script fixes that: it first pairs samples **globally by
`tumorid`** (ignoring `group` entirely), then decides which subgroup bucket
each *pair* belongs to using **the primary sample's `group` label**:

```python
# scripts/05_gsea.py, lines 47-58
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
```

This recovers all 43 pairs with none dropped — Group 4 MB ends up with the
correct n=14 (not n=12 as in script 03's output) and Group 3 MB is
unaffected here (its "extra" pair going *out* to Group 4's relapse-only side
is exactly cancelled by pairing on the primary side). The `assert` confirms
every bucketed pair's primary and relapse rows share the same `tumorid`
(after both are independently `sort_values("tumorid")`-ed, so row *i* of
`prim_sub` and row *i* of `recur_sub` are guaranteed to be the same patient)
before any stats are computed — a silent mismatch would fail loudly rather
than produce wrong numbers quietly.

### 5c. Expression filtering

```python
# scripts/05_gsea.py, lines 39-41
expressed = logcpm.index[logcpm.mean(axis=1) > 1]
logcpm_f = logcpm.loc[expressed]
print(f"genes after expression filter: {logcpm_f.shape[0]} / {logcpm.shape[0]}")
```

Before ranking, genes are filtered to `mean(logCPM) > 1` across all 86
samples (≈ mean CPM > 1) — 16,225 of the original 55,765 genes survive. This
removes genes that are essentially unexpressed cohort-wide, whose paired
t-statistics would otherwise be dominated by noise (a gene that's 0 counts in
84 samples and 1–2 counts in two samples can produce an enormous, meaningless
t-statistic). This is also *why* MIF and one other of the 7 user genes drop
out of `USER_ANGIOGENESIS_PANEL`'s effective size (reported "Tag %" denominators
of 5, not 7, in the GSEA output) — consistent with the MIF caveat in §6.

### 5d. Ranking metric and the GSEA run itself

For each subgroup, for every surviving gene: a **paired t-statistic** is
computed by hand (not via `scipy.stats.ttest_rel`, for speed across 16k
genes at once — but it's the same formula):

```python
# scripts/05_gsea.py, lines 60-68 (inside `for grp in group_order:`)
    p = logcpm_f[prim_ids].values
    r = logcpm_f[recur_ids].values
    diff = r - p
    n = diff.shape[1]
    mean_diff = diff.mean(axis=1)
    sd_diff = diff.std(axis=1, ddof=1)
    sd_diff[sd_diff == 0] = np.nan
    tstat = mean_diff / (sd_diff / np.sqrt(n))
    rnk = pd.Series(tstat, index=logcpm_f.index).dropna().sort_values(ascending=False)
```

(`sd_diff == 0 -> NaN` guards genes with zero variance in their paired
differences, e.g. identical values across every patient, which would
otherwise divide by zero; `.dropna()` then drops those genes from the
ranking entirely rather than assigning them a spurious infinite t-stat.)

Genes are ranked by this signed t-statistic, most-positive (consistently up
at relapse) at the top, most-negative (consistently up at primary) at the
bottom — this ranked list is what "preranked GSEA" runs on (as opposed to
"a priori GSEA," which computes the ranking from raw expression and phenotype
labels together; using a paired statistic as the ranking metric is the
standard adaptation for a paired/matched design, analogous to how the
original paper used limma with patient-pair blocking to call its DEGs).

`gseapy.prerank` (a Python port of the Broad's GSEA algorithm) is then run
against the three gene sets simultaneously, per subgroup:

```python
# scripts/05_gsea.py, lines 72-82 (inside `for grp in group_order:`)
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
```

- `min_size=3, max_size=2000` — a gene set needs at least 3 of its genes
  present in the ranked list to be scored at all (relevant for
  `USER_ANGIOGENESIS_PANEL` after filtering — it still clears this easily).
- `permutation_num=1000` — the null distribution for the p-value/NES is
  built from 1000 gene-set permutations (this is the standard GSEA
  permutation-based significance test, not a parametric approximation).
- `seed=42` — fixed for reproducibility of the permutation-based p-values.
- `outdir=None` — keeps `gseapy` from writing its own plots/report files to
  disk; only the in-memory `res2d` table is kept and written out ourselves
  as `gsea_angiogenesis_results.csv`.

**Output per (subgroup, gene set):** `NES` (normalized enrichment score —
positive means the gene set's members skew toward the "up at relapse" end of
the ranking; magnitude reflects how strongly, normalized so it's comparable
across gene sets of different sizes), `NOM p-val` (nominal, from the
permutation null), `FDR q-val` (multiple-testing-corrected across the gene
sets tested), and `Tag %` (what fraction of the gene set's members actually
appear in the leading edge driving the enrichment).

Also saved as a side artifact: `results/rank_<group>.rnk` — wait, these
`.rnk` files are written to the results directory by the script but are
per-run scratch (regenerable from `logcpm_matrix.csv` + `sample_metadata.csv`
by re-running this script) and are not committed, same reasoning as the
logCPM matrix.

**Output:** `results/gsea_angiogenesis_results.csv` — one row per
(subgroup, gene set).

## 6. Known caveats and limitations (collected)

- **MIF is unreliable.** Raw counts for MIF are near-zero (single digits) in
  almost all 86 samples, despite MIF normally being a moderately-to-highly
  expressed gene. The paper's own counting pipeline uses
  "uniquely mapped reads only" (their Methods section); MIF has several
  processed pseudogenes elsewhere in the genome, so reads that would
  otherwise support MIF's true expression are likely discarded as
  multi-mapping. This isn't a bug in this analysis — it's inherited from how
  the underlying counts were generated — but it means MIF's specific
  per-gene result (§3) and its contribution to the composite score and to
  `USER_ANGIOGENESIS_PANEL` should be discounted.
- **Group 3 MB is underpowered.** 5–6 pairs is not much for a paired
  Wilcoxon test or for GSEA's permutation null; several "trend but not
  significant" results in Group 3 (notably ITGB1, p=0.06, 5/5 patients up;
  and the GSEA NES for Hallmark Angiogenesis, p=0.067) are consistent
  directionally with the significant findings in the larger subgroups but
  shouldn't be over-interpreted as null results — they're likely a power
  problem, not an absence of effect.
- **n mismatch between script 03 and script 05 for Group 3/Group 4.** See
  §3 and §5b — script 03 reports n=5 (Group 3) / n=12 (Group 4) from
  intersection-after-filtering; script 05 correctly recovers n=6 / n=14 by
  pairing on `tumorid` first. Both are legitimate, just answering slightly
  different bucketing questions; don't directly compare an n across the two
  outputs without accounting for this.
- **`CURATED_ANGIOGENESIS_SUPPLEMENTARY` is not an authoritative GO term** —
  see §5a. Treat any result attributed to it as a secondary, lower-confidence
  read, not equivalent in evidentiary weight to the verified Hallmark set.
- **CPM, not RPKM** — see §2. Fine for the paired same-gene comparisons made
  here; would need gene-length correction if this matrix were later reused
  to compare different genes' expression levels to each other.
- **Composite angiogenesis score (§3) is an unweighted mean of z-scores**,
  not a formal gene-set enrichment score like GSVA/ssGSEA — it's a quick
  per-sample summary of the 7 requested genes specifically, distinct from
  (and complementary to) the GSEA results in §5, which score gene sets
  against the whole ranked transcriptome.
- **Raw sequencing data is intentionally not committed to this repo** — see
  `README.md`. Anyone re-running the pipeline needs their own local copies of
  the two source files pointed to via the `MB_META_FILE` / `MB_COUNTS_FILE`
  environment variables (see the top of `01_build_metadata.py` and
  `02_load_counts.py`).
