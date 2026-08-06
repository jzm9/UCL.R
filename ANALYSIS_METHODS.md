# scRNA-seq Analysis Methods & Justification

**Dataset:** GSE264326 — mouse leptomeningeal metastases (LPT_MET, n=4), primary
tumours (LPT_MET-matched, n=4), and non-tumour leptomeninges (LPT_WT, n=3).
11 samples total, 4 patient-matched MET/Primary pairs plus an unpaired
non-tumour baseline.

**Scripts referenced:** `R/scrna_leptomeningeal_analysis.R` (main Myriad
pipeline), `R/downstream_analysis/*.R` (local interrogation scripts).

This document walks through every analytical step taken, why it was chosen,
and — honestly, not defensively — where it deviates from best practice or
carries caveats that affect how the results should be interpreted. The
cluster identification section is treated in the most depth, since that step
carries the most subjective judgement and the most downstream consequence.

---

## 1. Data loading and QC

**What was done:** Each of the 11 samples loaded individually via `ReadMtx`
from GSM-specific matrix/features/barcodes files, `CreateSeuratObject` with
`min.cells = 3, min.features = 200`, then filtered per-sample on
`nFeature_RNA` (200–6000), `nCount_RNA` (>500), and `pct_mt` (<20%).

```r
samples <- data.frame(
  gsm        = c("GSM8216129", "GSM8216130", "GSM8216131",
                 "GSM8216132", "GSM8216133", "GSM8216134",
                 "GSM8216135", "GSM8216136", "GSM8216137",
                 "GSM8216138", "GSM8216139"),
  label      = c("LPT_MET_1", "Primary_1", "LPT_WT_1",
                 "LPT_MET_2", "Primary_2", "LPT_WT_2",
                 "LPT_MET_3", "Primary_3", "LPT_WT_3",
                 "LPT_MET_4", "Primary_4"),
  condition  = c("LPT_MET", "Primary", "LPT_WT",
                 "LPT_MET", "Primary", "LPT_WT",
                 "LPT_MET", "Primary", "LPT_WT",
                 "LPT_MET", "Primary"),
  patient    = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4),
  stringsAsFactors = FALSE
)

seurat_list <- lapply(seq_len(nrow(samples)), function(i) {
  s   <- samples[i, ]
  mtx <- ReadMtx(
    mtx      = file.path(data_dir, paste0(s$gsm, "_", s$label, "_matrix.mtx")),
    features = file.path(data_dir, paste0(s$gsm, "_", s$label, "_features.tsv")),
    cells    = file.path(data_dir, paste0(s$gsm, "_", s$label, "_barcodes.tsv"))
  )
  obj <- CreateSeuratObject(counts = mtx, project = s$label, min.cells = 3, min.features = 200)
  obj$sample    <- s$label
  obj$condition <- s$condition
  obj$patient   <- s$patient
  obj
})
names(seurat_list) <- samples$label

seurat_list <- lapply(seurat_list, function(obj) {
  obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  obj
})

seurat_list <- lapply(seurat_list, function(obj) {
  subset(obj,
    subset = nFeature_RNA > 200 &
             nFeature_RNA < 6000 &
             nCount_RNA   > 500  &
             pct_mt       < 20
  )
})
```

**Justification:** Per-sample filtering (rather than filtering after
merging) is correct practice — QC thresholds should account for
sample-to-sample technical variation, and filtering pre-merge avoids one
noisy sample dragging cells out of a shared threshold. The specific cutoffs
(200–6000 genes, <20% mito) are standard, permissive defaults for mouse
brain-adjacent tissue, not tuned to this dataset's own QC-metric
distributions.

**Caveat / not best practice:** The thresholds were not derived from this
dataset's own QC distributions (e.g. no per-sample MAD-based or knee-point
threshold, no visual inspection loop before committing to fixed cutoffs).
`qc_violin.pdf` is produced but the pipeline does not stop to let you adjust
thresholds per-sample based on it — it's descriptive, not adaptive. For a
paper-quality analysis, thresholds should be justified against the observed
distribution per sample (e.g. via `scater::isOutlier` or manual MAD-based
cutoffs), especially given 11 samples processed across what may be different
capture batches.

---

## 2. Normalisation (SCTransform)

**What was done:** `SCTransform(vars.to.regress = "pct_mt", variable.features.n = 3000)`
run **per sample**, before merging. Samples downsampled to a max of 5000
cells before normalisation, for memory reasons on the HPC node.

```r
max_cells <- 5000
seurat_list <- lapply(seurat_list, function(obj) {
  if (ncol(obj) > max_cells) obj <- obj[, sample(colnames(obj), max_cells)]
  obj <- SCTransform(obj, vars.to.regress = "pct_mt",
                     variable.features.n = 3000, verbose = FALSE)
  obj
})
```

**Justification:** SCTransform's regularized negative binomial regression is
current best practice over log-normalisation + scaling for UMI-based
scRNA-seq, and per-sample fitting is required before Harmony-style
integration, since SCTransform's regression should not be fit across
batches with different technical characteristics. Regressing out `pct_mt` is
a reasonable choice to reduce a known confound from cell stress/dying cells,
though it is one of several possible covariates.

**Caveat:** Cell-cycle score (`CellCycleScoring`) was **not** regressed out.
Given clusters 3/4/6 (in the all-cells object; see §4 below) turned out to
be pure cell-cycle signatures with no independent lineage marker in their
own top-10 differentially expressed genes, this omission has a direct,
visible consequence: proliferating cells clustered by cycle phase rather
than by cell type/lineage, which is exactly the artifact cell-cycle
regression exists to prevent. Whether to regress cell-cycle out is a
judgement call (it can also remove real biological signal, e.g. if
proliferation rate itself differs meaningfully between conditions — which is
plausible here, MET vs Primary), but the decision was not made deliberately;
it was simply not addressed. This should be treated as the single most
consequential unaddressed methodological choice in the whole pipeline.

The 5000-cell-per-sample downsampling cap is a pragmatic HPC memory
compromise, not a scientific choice — it discards real cells and could bias
against detecting rare populations if a sample was deep enough to exceed the
cap in a rare cell type. Not disclosed/quantified anywhere in the outputs
(no log of how many cells were dropped per sample).

---

## 3. Integration (Harmony) and clustering

**What was done:** PCA (50 PCs) on the SCT assay across the merged object,
`RunHarmony(group.by.vars = "sample", max_iter = 20)`, then UMAP/
`FindNeighbors`/`FindClusters` on the harmony-corrected embedding using 30
PCs and `resolution = 0.5` for the all-cells object (0.4/20 PCs for the
tumour-only re-clustering).

```r
merged <- merge(seurat_list[[1]], y = seurat_list[-1],
                add.cell.ids = names(seurat_list))
rm(seurat_list); gc()

DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- rownames(merged[["SCT"]]@scale.data)
merged <- RunPCA(merged, npcs = 50)

merged <- RunHarmony(merged, group.by.vars = "sample", max_iter = 20)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)
```

The tumour-only object is re-normalised and re-integrated from scratch after
subsetting (§4/§7), using the same pattern with `dims = 1:20` and
`resolution = 0.4`:

```r
tumour <- SCTransform(tumour, vars.to.regress = "pct_mt",
                      variable.features.n = 3000, verbose = FALSE)
DefaultAssay(tumour) <- "SCT"
VariableFeatures(tumour) <- rownames(tumour[["SCT"]]@scale.data)
tumour <- RunPCA(tumour, npcs = 50)
tumour <- RunHarmony(tumour, group.by.vars = "sample", max_iter = 20)
tumour <- RunUMAP(tumour, reduction = "harmony", dims = 1:20)
tumour <- FindNeighbors(tumour, reduction = "harmony", dims = 1:20)
tumour <- FindClusters(tumour, resolution = 0.4)
```

**Justification:** Harmony is a well-validated, commonly-used integration
method for this scale of dataset (11 samples, tens of thousands of cells)
and is appropriate given the goal is cross-sample comparison rather than
label transfer. Batch-correcting on `sample` (not `condition` or `patient`)
is correct — you want to remove technical/sample-of-origin variation while
preserving the biological (condition) variation you intend to test.

**Caveat / not best practice:** The clustering resolution (0.5 for all
cells, 0.4 for tumour-only) was **not derived by any systematic method** —
no `clustree` stability analysis, no silhouette/other resolution-selection
metric, just a fixed value. This was flagged explicitly earlier in this
project's development and the answer at the time was to skip clustree in
favour of direct marker-based cluster identification — a defensible
time-saving shortcut for exploratory analysis, but not something that should
be presented as a rigorously chosen resolution in a manuscript. The fact
that 23 clusters emerged, several of which collapse into essentially
redundant identities once labelled (see §4), is itself weak evidence the
resolution may be on the high side for the questions being asked, though it
also usefully over-splits things like the two endothelial subpopulations
that turned out to be biologically distinct (§4.3).

The number of PCs used (30 for all-cells, 20 for tumour-only) was not
validated against an elbow plot or `JackStraw`-type significance test in
this pipeline — it is a standard default, not a data-driven choice for this
dataset specifically.

---

## 4. Cluster identification — detailed treatment

This is the step given the most scrutiny, because it is the most subjective
and the most consequential for every downstream biological claim (TME
composition, "which cell types expand in MET").

### 4.1 Method used

1. `FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)`
   run on the all-cells, Harmony-integrated object (23 clusters at
   resolution 0.5), using the SCT assay (`PrepSCTFindMarkers()` run first
   to harmonise per-sample SCT models — required in Seurat v5 when multiple
   samples each carry their own SCT model, otherwise `FindMarkers`/
   `FindAllMarkers` throws an "unequal library sizes" error).

   ```r
   Idents(merged) <- "seurat_clusters"
   merged <- PrepSCTFindMarkers(merged)
   all_markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   write.csv(all_markers, "all_cells_cluster_markers_full.csv", row.names = FALSE)

   top10 <- all_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
   write.csv(top10, "all_cells_cluster_markers_top10.csv", row.names = FALSE)
   ```

2. Top-10 markers per cluster (by `avg_log2FC`) inspected manually.
3. Identities assigned by matching marker genes to literature/canonical
   cell-type markers, by eye, cluster by cluster, then applied to the object
   as a lookup table (`R/downstream_analysis/label_clusters.R`):

   ```r
   cluster_labels <- c(
     "0"  = "Tumour",
     "1"  = "Tumour",
     "2"  = "Smooth muscle",
     "3"  = "Tumour (S/G2M-phase)",
     "4"  = "Tumour (S/G2M-phase)",
     "5"  = "Pericytes",
     "6"  = "Tumour (S/G2M-phase)",
     "7"  = "Oligodendrocyte-like/glia",
     "8"  = "BBB endothelial",
     "9"  = "Schwann cells",
     "10" = "Fibroblasts",
     "11" = "Inflammatory macrophages/monocytes",
     "12" = "Microglia/macrophages",
     "13" = "Schwann cells",
     "14" = "Sex-linked/technical",
     "15" = "Tumour (stress/IEG state)",
     "16" = "Astrocytes",
     "17" = "Neutrophils",
     "18" = "T cells / NK cells",
     "19" = "Endothelial (Pecam1/Sox17-high)",
     "20" = "Homeostatic microglia",
     "21" = "Endothelial (Kdr/Mmrn2-high)",
     "22" = "Low-quality"
   )

   # unname() is required here: cluster_labels[as.character(...)] returns a
   # *named* vector (names = cluster IDs), and Seurat's `$<-` interprets a
   # named vector as needing name-based matching against cell barcodes,
   # which fails with "No cell overlap between new meta data and Seurat
   # object" if left named. unname() restores plain positional assignment.
   merged$cell_type <- factor(unname(cluster_labels[as.character(merged$seurat_clusters)]),
                               levels = unique(cluster_labels))

   merged$compartment <- ifelse(grepl("^Tumour", merged$cell_type), "Tumour", "TME")
   ```

4. Adjacent/overlapping identities (e.g. clusters 19 and 21, both initially
   labelled "Endothelial (vascular)") were revisited once their actual
   top-10 marker sets were checked: 19 is defined by `Pecam1`/`Sox17`
   (canonical pan-endothelial), 21 by `Kdr`/`Mmrn2`/`Bvht`/`Ch25h` — non-
   overlapping gene sets, split into separate labels rather than merged
   (reflected in the `cluster_labels` table above), since collapsing them
   would have hidden a real subpopulation distinction.

### 4.2 Is this best practice?

**Partially.** Manual, marker-gene-driven annotation against canonical
literature markers is a legitimate and widely-used approach — it is what
most published scRNA-seq papers actually do for the "final" annotation
layer, especially for well-characterised tissue compartments (immune,
vascular, stromal). The markers used here are, in the cases checked, genuine
canonical markers for their assigned identity (e.g. `Kcnj8`/`Abcc9`/`Rgs5`
for pericytes, `Cd3d/e/g` for T cells, `P2ry12`/`Tmem119`/`Cx3cr1` for
homeostatic microglia — all textbook-correct).

**Where it falls short of best practice:**

- **No automated cross-validation.** No reference-based label transfer
  (e.g. `SingleR`, Seurat's own label-transfer against a published
  reference atlas such as the Allen Brain or Tabula Muris) was used to
  independently corroborate the manual calls. Manual annotation alone is
  vulnerable to confirmation bias — you tend to find what you expect to
  find in a top-10 list, especially when curating literature markers
  yourself rather than testing them blind.
- **The manual annotation table was drafted by the user, then implemented
  in code without independent verification against the real
  `FindAllMarkers` output until explicitly requested.** When that
  cross-check was finally done (this session), it surfaced two genuine
  problems (§4.3, §4.4) — meaning the "trust but verify" step should have
  been the default order of operations, not an afterthought.
- **No per-cluster statistical confidence reported.** Cluster identity is
  presented as a categorical label with no accompanying certainty measure
  (e.g. proportion of top markers matching canonical signature, or a
  score like a module/gene-set-enrichment p-value per cluster per
  candidate identity). A cluster like #14 (see below) should arguably not
  have been given a confident-sounding label at all.
- **No doublet detection step** (e.g. `DoubletFinder`, `scDblFinder`)
  anywhere in the pipeline before clustering. In a dataset combining tumour
  and non-tumour tissue with many rare TME populations, doublets between
  tumour cells and TME cells are a real risk and could plausibly explain
  some of the CopyKAT tumour/TME misclassification described in §5.

### 4.3 Specific findings from the cross-check (this session)

Cross-checking the stated cluster→marker table against the actual
`all_cells_cluster_markers_top10.csv` output found:

| Cluster(s) | Assigned identity | Verdict | Note |
|---|---|---|---|
| 3, 4, 6 | "Tumour (S/G2M-phase)" | **Weak support** | Real top-10 markers are *purely* mitotic/S-phase genes (`Ccnb1/2`, `Cenpf`, `Cdc20` for 3; replication histones for 4; `Mcm2-6` for 6). None show tumour-lineage markers (`Neurod1`/`Atoh1`/`Barhl1`) in their own top-10. The tumour call rests on UMAP-proximity inference to clusters 0/1, not on independent marker evidence — this is precisely the ambiguity that cell-cycle regression (§2) is designed to resolve, and its absence is the direct cause. |
| 13 | "Schwann cells" | **Partially supported, different genes** | Only `Scn7a` of the three markers originally cited (`Mag`/`Pllp`/`Ncmap`) actually appears in cluster 13's real top-10; those three are cluster 9's markers. Cluster 13's own top-10 does include `Cdh19`, a genuine independent Schwann-cell marker, so the identity call still holds — but the original justification was wrong even though the conclusion happened to be right. |
| 14 | "Sex-linked/technical" | **Questionable as a discrete identity** | `Xist`/`Tsix` are present but so is `Gria2` (glutamate receptor, neuronal) and `Kcnq1ot1`/`Miat` (imprinted lncRNAs unrelated to sex). This looks more like a low-specificity/mixed cluster than a genuine "technical" artifact category — treating it as a real, dismissable cell-type category (as opposed to a cluster requiring further splitting or QC review) is not well supported by the marker evidence. |
| 19 vs 21 | Originally both "Endothelial (vascular)"; **split** during this cross-check into "Endothelial (Pecam1/Sox17-high)" and "Endothelial (Kdr/Mmrn2-high)" | **Correction applied** | Merging these would have hidden a real subpopulation split; supported by non-overlapping top-10 marker sets. |
| 22 | "Low-quality, ignore" | **Supported** | Several markers have `p_val_adj = 1` (non-significant), consistent with a low-information/junk cluster. Appropriately excluded from biological interpretation. |

All other cluster identities (0, 1, 2, 5, 7, 8, 9, 10, 11, 12, 15, 16, 17,
18, 20) were checked and their stated canonical markers do appear, with
high specificity (low `pct.2`), in the real top-10 output — these are
well-supported.

### 4.4 The `cell_type` assignment bug (fixed, but worth flagging methodologically)

The first implementation of `label_clusters.R` assigned labels via:

```r
merged$cell_type <- factor(cluster_labels[as.character(merged$seurat_clusters)], ...)
```

This produced `Error: No cell overlap between new meta data and Seurat object`,
because `cluster_labels[as.character(...)]` returns a *named* vector (names
= cluster IDs, not cell barcodes), and Seurat's `$<-` assignment tries to
match those names against cell barcodes rather than assigning positionally.
Fixed by wrapping in `unname()`. This is purely a coding bug, not a
methodological one, but it's included here because it's exactly the kind of
silent-until-caught error that supports the case for automated
cross-validation (§4.2) — a bug that would have produced *no* cell-type
metadata at all if it had failed differently (e.g. partial matching instead
of a hard error).

---

## 5. Tumour/TME classification (CopyKAT)

**What was done:** CopyKAT run **per tumour sample** (LPT_MET/Primary),
each time using the pooled LPT_WT (non-tumour leptomeninges) cells as the
explicit normal reference (`norm.cell.names`), with `genome = "mm10"`
(required for mouse — CopyKAT defaults to human hg20 coordinates, and
without this argument every gene fails coordinate annotation, producing an
"all cells are filtered" error, encountered and fixed during development).
Cells called `aneuploid` and drawn from an LPT_MET/Primary sample are
labelled tumour.

```r
normal_barcodes <- WhichCells(merged, expression = condition == "LPT_WT")
tumour_samples  <- unique(merged$sample[merged$condition %in% c("LPT_MET", "Primary")])

copykat_preds <- lapply(tumour_samples, function(samp) {
  cells_in_sample <- WhichCells(merged, expression = sample == samp)
  cells_to_use    <- c(cells_in_sample, normal_barcodes)
  raw             <- as.matrix(GetAssayData(merged, layer = "counts")[, cells_to_use])

  out_dir <- file.path(data_dir, paste0("copykat_", samp))
  dir.create(out_dir, showWarnings = FALSE)
  old_wd <- getwd(); setwd(out_dir)

  ck <- copykat(
    rawmat          = raw,
    norm.cell.names = normal_barcodes,
    id.type         = "S",
    genome          = "mm10",   # required for mouse - defaults to human hg20 otherwise,
                                 # which causes an "all cells are filtered" error since no
                                 # mouse gene symbol maps to a chromosome coordinate
    ngene.chr       = 5,
    win.size        = 25,
    KS.cut          = 0.1,
    sam.name        = samp,
    distance        = "euclidean",
    n.cores         = 8
  )
  setwd(old_wd)
  ck$prediction
})

all_preds <- do.call(rbind, copykat_preds)
merged$copykat_pred <- NA
merged$copykat_pred[match(all_preds$cell.names, colnames(merged))] <- all_preds$copykat.pred

merged$is_tumour <- merged$copykat_pred == "aneuploid" &
                    merged$condition %in% c("LPT_MET", "Primary")
```

**Justification:** Using genuine non-tumour tissue from the same
anatomical/experimental context (LPT_WT) as the CopyKAT normal reference is
the correct approach — it's a much stronger reference than CopyKAT's default
of inferring "diploid-like" cells internally from the same sample, which is
less reliable when tumour purity is high.

**Caveat — the most significant known limitation of the whole pipeline:**
Cluster-marker-based inspection of the CopyKAT "tumour" call (performed
earlier in this project) found approximately **14% of CopyKAT-classified
"tumour" cells were actually TME cell types** (macrophages, microglia,
astrocytes, fibroblasts/pericytes, smooth muscle, dendritic cells, Schwann
cells, neutrophils — clusters 7–14 in the tumour-only re-clustering carried
these TME identities despite being CopyKAT-"tumour"-classified). This is a
real, non-trivial contamination rate. CNV-inference methods like CopyKAT are
known to underperform on tumours with low aneuploidy burden or on stromal
cells with subclonal/technical CNV-like noise; this is a known limitation of
the method generally, not unique to this implementation, but 14% is high
enough that:

- **Any DE result or angiogenesis-pathway conclusion drawn from the
  "tumour-only" object (`seurat_tumour_cells.rds`) should be treated
  cautiously** — several TME marker genes (e.g. `Trem1`, `Lrg1`, `Thbd`,
  `Col3a1`, `Kcnj8`, `Pglyrp1`) showed apparently "significant" differential
  expression in the tumour-only pseudobulk DE, which is very plausibly an
  artifact of this contamination (TME cell composition shifting between
  MET/Primary, mislabelled as "tumour" signal) rather than genuine
  tumour-intrinsic biology.
- This is why the pipeline was restructured mid-project to **also** save
  and separately analyse the full all-cells object
  (`seurat_all_cells_integrated.rds`) with `is_tumour` retained as
  metadata, rather than working exclusively from the CopyKAT-filtered
  subset — the TME-composition analysis (§6) is deliberately run on this
  unfiltered object for that reason.
- **Best-practice alternative not implemented:** a second, orthogonal
  tumour-calling method (e.g. `inferCNV`, or a marker-gene-score-based
  classifier using known tumour-intrinsic genes as a sanity check against
  CopyKAT) was not run in parallel. Cross-validating CNV-based calls against
  an independent method is standard practice when tumour purity/contamination
  is a known concern, and was not done here.

---

## 6. TME composition analysis

**What was done:** Run on the full all-cells object (not the CopyKAT-filtered
subset, for the reasons in §5). Cluster (later, labelled cell-type)
proportions computed per patient/condition; paired Wilcoxon signed-rank test
(MET vs Primary, paired by patient, n=4 pairs) used to test for proportional
shifts; LPT_WT included separately as an unpaired reference baseline, not
part of the formal paired test.

**Justification:** Paired testing is correct given the study design
(matched MET/Primary pairs per patient) — an unpaired test would discard
information and inflate variance unnecessarily. Excluding LPT_WT from the
paired test (rather than trying to force it into a three-group paired
design it doesn't fit) is the right call, since it has no patient-matched
partner.

```r
props <- merged@meta.data %>%
  filter(condition %in% c("LPT_MET", "Primary")) %>%
  count(patient, condition, cell_type) %>%
  group_by(patient, condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_wide <- props %>%
  select(patient, condition, cell_type, prop) %>%
  tidyr::pivot_wider(names_from = condition, values_from = prop, values_fill = 0)

cell_type_prop_test <- prop_wide %>%
  group_by(cell_type) %>%
  summarise(
    mean_MET     = mean(LPT_MET, na.rm = TRUE),
    mean_Primary = mean(Primary, na.rm = TRUE),
    p_value      = tryCatch(wilcox.test(LPT_MET, Primary, paired = TRUE)$p.value,
                             error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  arrange(p_value)
```

**Caveats:**

- **Compositional-data problem, acknowledged but not corrected for.**
  Cluster/cell-type proportions are computed within each
  patient/condition and necessarily sum to 1. This means an increase in
  one population's proportion mechanically forces others down, even absent
  any real change in that other population's abundance. No compositional
  data method (e.g. centered log-ratio transform, `scCODA`,
  `propeller`) was used to correct for this — proportions were compared
  directly with a standard Wilcoxon test on raw proportions. This is a
  common practice in the field but is increasingly recognised as suboptimal;
  a compositionally-aware method would be more defensible for a final
  manuscript-quality claim about which populations "genuinely" expand or
  contract.
- **Statistical floor at n=4 patient pairs.** With 4 matched pairs, the
  minimum achievable two-sided exact p-value from a paired Wilcoxon test is
  0.125 — meaning "not statistically significant" for most of the tested
  populations at this sample size does not equate to "no real effect." This
  was noted during the analysis but bears repeating here: any population
  showing a consistent, large-magnitude shift across all 4 pairs should be
  taken seriously as a plausible biological signal even without reaching a
  conventional significance threshold, given the test's structural
  ceiling.
- **No multiple-testing correction applied** across the ~14-23 populations
  tested simultaneously for proportional shift. At this many parallel
  tests, some nominal p-values under 0.05 should be expected by chance
  alone; FDR correction (e.g. Benjamini-Hochberg) was not applied to the
  `cluster_proportion_test_MET_vs_Primary.csv` /
  `cell_type_proportion_test_MET_vs_Primary.csv` outputs.

---

## 7. Differential expression (pseudobulk DESeq2)

**What was done:** `AggregateExpression()` used to sum raw RNA counts per
sample/condition/patient (pseudobulk), then `FindMarkers(test.use = "DESeq2")`
comparing LPT_MET vs Primary pseudo-samples (n=4 vs n=4, tumour cells only).

```r
DefaultAssay(tumour) <- "RNA"
tumour_met_prim <- subset(tumour, subset = condition %in% c("LPT_MET", "Primary"))

pseudobulk <- AggregateExpression(tumour_met_prim,
  group.by = c("sample", "condition", "patient"),
  assays    = "RNA",
  return.seurat = TRUE
)

Idents(pseudobulk) <- "condition"
# Seurat replaces underscores with dashes in identity class labels
# (e.g. "LPT_MET" -> "LPT-MET"), so look up the actual sanitised levels
# rather than hardcoding the original condition strings - hardcoding
# "LPT_MET" here previously failed with "Cannot find the following
# identities in the object: LPT_MET".
ident_levels  <- levels(Idents(pseudobulk))
ident_met     <- grep("MET", ident_levels, value = TRUE)
ident_primary <- grep("Primary", ident_levels, value = TRUE)

de_results <- FindMarkers(
  pseudobulk,
  ident.1    = ident_met,
  ident.2    = ident_primary,
  test.use   = "DESeq2",
  min.pct    = 0.1,
  logfc.threshold = 0
)
```

**Justification:** Pseudobulk aggregation before DESeq2 is correct and
current best practice for scRNA-seq DE testing between conditions —
testing directly on single-cell-level pseudo-replicates (treating each cell
as an independent observation) is a well-known statistical error that
massively inflates false-positive rates by ignoring within-sample
cell-to-cell correlation; pseudobulk avoids this by returning the unit of
replication to the biological sample, matching the true n=4 vs n=4 design.

**Caveats:**

- As noted in §5, this DE analysis was run on the CopyKAT-filtered
  tumour-only object, which carries the ~14% TME contamination — so several
  "significant" DE genes in this output are plausibly composition artifacts
  rather than tumour-intrinsic expression changes. This is the single
  biggest interpretive caveat on the DE results specifically.
- Only n=4 vs n=4 pseudo-samples — DESeq2's dispersion estimation is
  reasonably robust at this n but still limited; results should be treated
  as hypothesis-generating rather than confirmatory without independent
  validation (e.g. qPCR/IHC follow-up on top hits).

---

## 8. GSEA / pathway analysis

**What was done:** `fgsea` run against combined Hallmark + KEGG + Reactome
mouse gene sets from `msigdbr`, ranked by `avg_log2FC × -log10(p_val)` from
the pseudobulk DE. Angiogenesis-related gene sets specifically extracted via
keyword grep (`ANGIOGEN|VEGF|NOTCH|HIF|HYPOXIA|VESSEL|VASCULO`) for focused
follow-up (enrichment plots, marker-gene heatmap, module score).

**Justification:** `fgsea` is a standard, appropriately fast GSEA
implementation; combining Hallmark/KEGG/Reactome is reasonable pathway
coverage. Ranking by a combined effect-size/significance metric (rather
than log2FC alone) is a defensible choice that down-weights large but
noisy fold-changes.

```r
de_results$rank_metric <- de_results$avg_log2FC * -log10(de_results$p_val + 1e-300)
ranked <- setNames(de_results$rank_metric, de_results$gene)
ranked <- sort(ranked, decreasing = TRUE)

hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)
all_genesets <- c(hallmark, kegg, reactome)

set.seed(42)
gsea_res <- fgsea(pathways = all_genesets, stats = ranked,
                   minSize = 10, maxSize = 500, nPermSimple = 10000)

# Focused angiogenesis follow-up: keyword filter over the unbiased GSEA result
angio_sets <- grep("ANGIOGEN|VEGF|NOTCH|HIF|HYPOXIA|VESSEL|VASCULO",
                   names(all_genesets), value = TRUE, ignore.case = TRUE)
angio_gsea <- gsea_res[gsea_res$pathway %in% angio_sets, ]

core_angio_genes <- unique(unlist(all_genesets[angio_sets]))
core_angio_genes <- core_angio_genes[core_angio_genes %in% rownames(tumour)]
avg_expr <- AverageExpression(tumour_met_prim, features = core_angio_genes,
                               group.by = "sample", layer = "data")$RNA
avg_z <- t(scale(t(avg_expr)))

hallmark_angio <- all_genesets[["HALLMARK_ANGIOGENESIS"]]
hallmark_angio <- hallmark_angio[hallmark_angio %in% rownames(tumour)]
tumour <- AddModuleScore(tumour, features = list(hallmark_angio),
                          name = "angio_score", ctrl = 100)
```

**Caveats:**

- The angiogenesis gene-set selection via keyword grep is a convenience
  filter, not a pre-registered hypothesis — it's appropriate for
  hypothesis-generating exploration (which is how it's used here) but
  should not be presented as if angiogenesis was the sole a priori focus of
  an unbiased genome-wide GSEA; the GSEA itself is unbiased, but the
  reporting narrows to angiogenesis afterward.
- GSEA results inherit the same tumour/TME-contamination caveat as the
  underlying DE ranks (§5, §7).
- The angiogenesis module score (`AddModuleScore`) uses only the Hallmark
  `HALLMARK_ANGIOGENESIS` gene set, not the union of all angiogenesis-related
  sets identified by the broader keyword grep — a scope inconsistency
  between the heatmap (uses the broader keyword-matched gene set) and the
  module score (uses only Hallmark) that isn't flagged anywhere in the
  script's comments.

---

## 9. Overall summary of best-practice gaps

Ranked roughly by how much each affects the reliability of the biological
conclusions:

1. **CopyKAT tumour/TME contamination (~14%)** flowing into the tumour-only
   DE and GSEA results — the single largest source of potential
   misinterpretation. Mitigated (not eliminated) by also running TME
   composition analysis on the unfiltered object.
2. **No cell-cycle regression**, directly causing clusters 3/4/6 to be
   defined by cell-cycle phase rather than lineage, with the tumour
   identity for those clusters resting on UMAP-proximity inference rather
   than independent marker evidence.
3. **Manual cluster annotation with no independent/automated
   cross-validation** (e.g. reference-based label transfer), and no
   doublet detection step before clustering.
4. **No compositional-data-aware statistics** for the TME proportion
   analysis, and no multiple-testing correction across the many populations
   tested.
5. **Clustering resolution and PC count chosen by convention, not
   data-driven selection** (no clustree, no elbow/JackStraw validation).
6. **QC thresholds fixed rather than derived from this dataset's own
   distributions.**

None of these gaps invalidate the pipeline — the core methodological choices
(SCTransform, Harmony, pseudobulk DESeq2, fgsea) are all appropriate and
current best practice at the architecture level. The gaps are in the details
of parameter/threshold selection and in validation steps that were skipped
for development speed during an iterative, exploratory HPC debugging
process. For any result destined for publication, items 1–4 above should be
addressed before the findings are treated as confirmatory rather than
exploratory.

---

## Appendix A: Full pipeline script (`R/scrna_leptomeningeal_analysis.R`)

Complete, unabridged, current version — every line as run on Myriad.

```r
# scRNA-seq Analysis: Leptomeningeal Metastases vs Primary Tumours
# GSE264326 — 11 samples: LPT_MET (n=4), Primary (n=4), LPT_WT (n=3)
# Focus: tumour cells only, angiogenesis pathways, LARGE1, unbiased GSEA

# ── 0. Install / load packages ────────────────────────────────────────────────
packages <- c(
  "Seurat", "harmony", "ggplot2", "dplyr", "patchwork",
  "msigdbr", "DESeq2",
  "fgsea", "ggrepel", "viridis", "ComplexHeatmap",
  "Matrix", "BiocParallel", "glmGamPoi"
)
installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install)) {
  if (!"BiocManager" %in% installed) install.packages("BiocManager")
  BiocManager::install(to_install, ask = FALSE)
}
invisible(lapply(packages, library, character.only = TRUE))

# CopyKAT: install from GitHub if needed (no JAGS dependency)
if (!"copykat" %in% rownames(installed.packages())) {
  if (!"devtools" %in% rownames(installed.packages())) install.packages("devtools")
  devtools::install_github("navinlabcode/copykat")
}
library(copykat)

# ── 1. Paths & sample metadata ───────────────────────────────────────────────
data_dir <- "/home/sejkmor/Scratch/taylordataset/taylor25"

samples <- data.frame(
  gsm        = c("GSM8216129", "GSM8216130", "GSM8216131",
                 "GSM8216132", "GSM8216133", "GSM8216134",
                 "GSM8216135", "GSM8216136", "GSM8216137",
                 "GSM8216138", "GSM8216139"),
  label      = c("LPT_MET_1", "Primary_1", "LPT_WT_1",
                 "LPT_MET_2", "Primary_2", "LPT_WT_2",
                 "LPT_MET_3", "Primary_3", "LPT_WT_3",
                 "LPT_MET_4", "Primary_4"),
  condition  = c("LPT_MET", "Primary", "LPT_WT",
                 "LPT_MET", "Primary", "LPT_WT",
                 "LPT_MET", "Primary", "LPT_WT",
                 "LPT_MET", "Primary"),
  patient    = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4),
  stringsAsFactors = FALSE
)

# ── 2. Load raw counts ────────────────────────────────────────────────────────
seurat_list <- lapply(seq_len(nrow(samples)), function(i) {
  s   <- samples[i, ]
  mtx <- ReadMtx(
    mtx      = file.path(data_dir, paste0(s$gsm, "_", s$label, "_matrix.mtx")),
    features = file.path(data_dir, paste0(s$gsm, "_", s$label, "_features.tsv")),
    cells    = file.path(data_dir, paste0(s$gsm, "_", s$label, "_barcodes.tsv"))
  )
  obj <- CreateSeuratObject(counts = mtx, project = s$label, min.cells = 3, min.features = 200)
  obj$sample    <- s$label
  obj$condition <- s$condition
  obj$patient   <- s$patient
  obj
})
names(seurat_list) <- samples$label

# ── 3. QC & filtering ─────────────────────────────────────────────────────────
seurat_list <- lapply(seurat_list, function(obj) {
  obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  obj
})

# QC summary table (no merge needed — avoids RAM spike)
qc_summary <- do.call(rbind, lapply(names(seurat_list), function(nm) {
  obj <- seurat_list[[nm]]
  data.frame(
    sample        = nm,
    n_cells_raw   = ncol(obj),
    median_genes  = median(obj$nFeature_RNA),
    median_counts = median(obj$nCount_RNA),
    median_pct_mt = median(obj$pct_mt)
  )
}))
write.csv(qc_summary, file.path(data_dir, "qc_summary_pre_filter.csv"), row.names = FALSE)

# QC violin plots per sample individually (low memory)
pdf(file.path(data_dir, "qc_violin.pdf"), width = 8, height = 10)
for (nm in names(seurat_list)) {
  p <- VlnPlot(seurat_list[[nm]],
               features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
               pt.size = 0, ncol = 3) &
    ggtitle(nm) &
    theme(axis.text.x = element_blank())
  print(p)
}
dev.off()

# Filter per sample
seurat_list <- lapply(seurat_list, function(obj) {
  subset(obj,
    subset = nFeature_RNA > 200 &
             nFeature_RNA < 6000 &
             nCount_RNA   > 500  &
             pct_mt       < 20
  )
})
gc()

# ── 4. SCTransform per sample (normalise + variable features + regression) ────
# SCTransform is far more memory-efficient than ScaleData with vars.to.regress.
# Downsample to max 5000 cells per sample to keep RAM manageable.
max_cells <- 5000
seurat_list <- lapply(seurat_list, function(obj) {
  if (ncol(obj) > max_cells) obj <- obj[, sample(colnames(obj), max_cells)]
  obj <- SCTransform(obj, vars.to.regress = "pct_mt",
                     variable.features.n = 3000, verbose = FALSE)
  obj
})

# ── 5. Merge & integrate with Harmony ────────────────────────────────────────
merged <- merge(seurat_list[[1]], y = seurat_list[-1],
                add.cell.ids = names(seurat_list))
rm(seurat_list); gc()

# Use SCT assay for PCA; PrepSCTFindMarkers harmonises residuals across samples
DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- rownames(merged[["SCT"]]@scale.data)
merged <- RunPCA(merged, npcs = 50)

merged <- RunHarmony(merged, group.by.vars = "sample", max_iter = 20)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)

pdf(file.path(data_dir, "umap_all_cells.pdf"), width = 14, height = 6)
p1 <- DimPlot(merged, group.by = "condition", label = FALSE) + ggtitle("Condition")
p2 <- DimPlot(merged, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
print(p1 | p2)
dev.off()

# ── 6. Tumour cell identification via CopyKAT ─────────────────────────────────
# CopyKAT uses a Bayesian model to infer CNV from scRNA counts.
# LPT_WT cells (non-tumour leptomeningeal) serve as the normal reference.
# Run CopyKAT per tumour sample, pooling LPT_WT normals as reference each time.

normal_barcodes <- WhichCells(merged, expression = condition == "LPT_WT")
tumour_samples  <- unique(merged$sample[merged$condition %in% c("LPT_MET", "Primary")])

copykat_preds <- lapply(tumour_samples, function(samp) {
  cells_in_sample <- WhichCells(merged, expression = sample == samp)
  cells_to_use    <- c(cells_in_sample, normal_barcodes)
  raw             <- as.matrix(GetAssayData(merged, layer = "counts")[, cells_to_use])

  out_dir <- file.path(data_dir, paste0("copykat_", samp))
  dir.create(out_dir, showWarnings = FALSE)
  old_wd <- getwd(); setwd(out_dir)

  ck <- copykat(
    rawmat          = raw,
    norm.cell.names = normal_barcodes,
    id.type         = "S",
    genome          = "mm10",
    ngene.chr       = 5,
    win.size        = 25,
    KS.cut          = 0.1,
    sam.name        = samp,
    distance        = "euclidean",
    n.cores         = 8
  )
  setwd(old_wd)
  ck$prediction
})

all_preds <- do.call(rbind, copykat_preds)
merged$copykat_pred <- NA
merged$copykat_pred[match(all_preds$cell.names, colnames(merged))] <- all_preds$copykat.pred

# Tumour cells = aneuploid in CopyKAT, from LPT_MET or Primary samples
merged$is_tumour <- merged$copykat_pred == "aneuploid" &
                    merged$condition %in% c("LPT_MET", "Primary")
merged$cnv_score <- as.numeric(merged$copykat_pred == "aneuploid")

pdf(file.path(data_dir, "umap_tumour_classification.pdf"), width = 12, height = 5)
p1 <- FeaturePlot(merged, features = "cnv_score", cols = c("grey90", "red3")) +
  ggtitle("Inferred CNV score")
p2 <- DimPlot(merged, cells.highlight = WhichCells(merged, expression = is_tumour == TRUE),
              cols.highlight = "red3", cols = "grey85") +
  ggtitle("Classified tumour cells")
print(p1 | p2)
dev.off()

# ── 7. Subset to tumour cells ─────────────────────────────────────────────────
tumour <- subset(merged, subset = is_tumour == TRUE)

# Save the full all-cells object here, before freeing it - step 13 previously
# tried to saveRDS(merged, ...) after rm(merged), which fails with
# "object 'merged' not found".
saveRDS(merged, file.path(data_dir, "seurat_all_cells_integrated.rds"))
rm(merged); gc()  # free the full merged object now that tumour subset is done

tumour <- SCTransform(tumour, vars.to.regress = "pct_mt",
                      variable.features.n = 3000, verbose = FALSE)
DefaultAssay(tumour) <- "SCT"
VariableFeatures(tumour) <- rownames(tumour[["SCT"]]@scale.data)
tumour <- RunPCA(tumour, npcs = 50)
tumour <- RunHarmony(tumour, group.by.vars = "sample", max_iter = 20)
tumour <- RunUMAP(tumour, reduction = "harmony", dims = 1:20)
tumour <- FindNeighbors(tumour, reduction = "harmony", dims = 1:20)
tumour <- FindClusters(tumour, resolution = 0.4)

pdf(file.path(data_dir, "umap_tumour_cells.pdf"), width = 14, height = 6)
p1 <- DimPlot(tumour, group.by = "condition", cols = c("LPT_MET" = "#E63946", "Primary" = "#457B9D")) +
  ggtitle("Tumour cells by condition")
p2 <- DimPlot(tumour, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Tumour clusters")
p3 <- DimPlot(tumour, group.by = "patient", label = FALSE) +
  ggtitle("Patient")
print((p1 | p2 | p3))
dev.off()

# ── 8. LARGE1 expression ──────────────────────────────────────────────────────
pdf(file.path(data_dir, "LARGE1_expression.pdf"), width = 14, height = 5)
p1 <- FeaturePlot(tumour, features = "Large1", cols = c("grey90", "darkblue")) +
  ggtitle("Large1 — UMAP")
p2 <- VlnPlot(tumour, features = "Large1", group.by = "condition",
              cols = c("LPT_MET" = "#E63946", "Primary" = "#457B9D"), pt.size = 0) +
  ggtitle("Large1 expression by condition")
p3 <- ggplot(tumour@meta.data, aes(x = condition, y = FetchData(tumour, vars = "Large1")[, 1], fill = condition)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = c("LPT_MET" = "#E63946", "Primary" = "#457B9D")) +
  labs(y = "Large1 (normalised)", x = NULL) +
  theme_classic() + theme(legend.position = "none")
print(p1 | p2 | p3)
dev.off()

# Large1 statistical test (Wilcoxon, paired by patient)
large1_expr <- FetchData(tumour, vars = c("Large1", "condition", "patient"))
large1_summary <- large1_expr %>%
  group_by(condition, patient) %>%
  summarise(mean_LARGE1 = mean(Large1), .groups = "drop")
large1_wide <- tidyr::pivot_wider(large1_summary, names_from = condition, values_from = mean_LARGE1)
large1_test  <- wilcox.test(large1_wide$LPT_MET, large1_wide$Primary, paired = TRUE)
message("LARGE1 Wilcoxon test (paired): p = ", round(large1_test$p.value, 4))
write.csv(large1_wide, file.path(data_dir, "LARGE1_per_patient.csv"), row.names = FALSE)

# ── 9. Pseudobulk DE: LPT_MET vs Primary ─────────────────────────────────────
# Aggregate counts per sample (pseudobulk) to avoid inflation of n
DefaultAssay(tumour) <- "RNA"
tumour_met_prim <- subset(tumour, subset = condition %in% c("LPT_MET", "Primary"))

pseudobulk <- AggregateExpression(tumour_met_prim,
  group.by = c("sample", "condition", "patient"),
  assays    = "RNA",
  return.seurat = TRUE
)

Idents(pseudobulk) <- "condition"
# Seurat replaces underscores with dashes in identity class labels
# (e.g. "LPT_MET" -> "LPT-MET"), so look up the actual sanitised levels
# rather than hardcoding the original condition strings.
ident_levels <- levels(Idents(pseudobulk))
ident_met     <- grep("MET", ident_levels, value = TRUE)
ident_primary <- grep("Primary", ident_levels, value = TRUE)

de_results <- FindMarkers(
  pseudobulk,
  ident.1    = ident_met,
  ident.2    = ident_primary,
  test.use   = "DESeq2",
  min.pct    = 0.1,
  logfc.threshold = 0
)
de_results$gene <- rownames(de_results)
de_results <- de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]
write.csv(de_results, file.path(data_dir, "DE_LPT_MET_vs_Primary.csv"), row.names = TRUE)

# Volcano plot
de_results$label <- ifelse(
  (abs(de_results$avg_log2FC) > 1 & de_results$p_val_adj < 0.05) |
    de_results$gene == "Large1",
  de_results$gene, ""
)
de_results$colour <- case_when(
  de_results$gene == "Large1"                                               ~ "LARGE1",
  de_results$avg_log2FC >  1 & de_results$p_val_adj < 0.05                 ~ "Up in MET",
  de_results$avg_log2FC < -1 & de_results$p_val_adj < 0.05                 ~ "Up in Primary",
  TRUE                                                                       ~ "NS"
)
pdf(file.path(data_dir, "volcano_MET_vs_Primary.pdf"), width = 8, height = 7)
ggplot(de_results, aes(avg_log2FC, -log10(p_val_adj), colour = colour, label = label)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text_repel(size = 3, max.overlaps = 30) +
  scale_colour_manual(values = c("Up in MET" = "#E63946", "Up in Primary" = "#457B9D",
                                 "LARGE1" = "darkgreen", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  labs(x = "log2 FC (MET / Primary)", y = "-log10 adj. p-value",
       title = "LPT_MET vs Primary — tumour cells (pseudobulk DESeq2)") +
  theme_classic()
dev.off()

# ── 10. Gene set enrichment (GSEA) ────────────────────────────────────────────
# Ranked gene list by log2FC × -log10(p)
de_results$rank_metric <- de_results$avg_log2FC * -log10(de_results$p_val + 1e-300)
ranked <- setNames(de_results$rank_metric, de_results$gene)
ranked <- sort(ranked, decreasing = TRUE)

# MSigDB gene sets: Hallmark + KEGG
hallmark <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

kegg <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

all_genesets <- c(hallmark, kegg, reactome)

set.seed(42)
gsea_res <- fgsea(
  pathways  = all_genesets,
  stats     = ranked,
  minSize   = 10,
  maxSize   = 500,
  nPermSimple = 10000
)
gsea_res <- gsea_res[order(gsea_res$NES, decreasing = TRUE), ]
write.csv(gsea_res[, -8], file.path(data_dir, "GSEA_MET_vs_Primary.csv"), row.names = FALSE)  # col 8 = leadingEdge list

# Top pathways plot (top 20 up + top 20 down)
top_paths <- rbind(
  head(gsea_res[gsea_res$padj < 0.05 & gsea_res$NES > 0, ], 20),
  head(gsea_res[gsea_res$padj < 0.05 & gsea_res$NES < 0, ], 20)
)
top_paths$pathway_short <- gsub("HALLMARK_|KEGG_|REACTOME_", "", top_paths$pathway)
top_paths$pathway_short <- stringr::str_trunc(top_paths$pathway_short, 55)
top_paths$direction <- ifelse(top_paths$NES > 0, "Up in MET", "Up in Primary")

pdf(file.path(data_dir, "GSEA_top_pathways.pdf"), width = 11, height = 10)
ggplot(top_paths, aes(x = NES, y = reorder(pathway_short, NES), fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up in MET" = "#E63946", "Up in Primary" = "#457B9D")) +
  geom_vline(xintercept = 0, colour = "black") +
  labs(x = "Normalised Enrichment Score", y = NULL,
       title = "GSEA: LPT_MET vs Primary (tumour cells, Hallmark + KEGG + Reactome)") +
  theme_classic() + theme(legend.title = element_blank())
dev.off()

# ── 11. Angiogenesis pathways — deep dive ────────────────────────────────────
angio_sets <- grep("ANGIOGEN|VEGF|NOTCH|HIF|HYPOXIA|VESSEL|VASCULO",
                   names(all_genesets), value = TRUE, ignore.case = TRUE)
message("Angiogenesis-related gene sets found: ", paste(angio_sets, collapse = ", "))

angio_gsea <- gsea_res[gsea_res$pathway %in% angio_sets, ]
write.csv(angio_gsea[, -8], file.path(data_dir, "GSEA_angiogenesis_pathways.csv"), row.names = FALSE)

# Enrichment plots for significant angiogenesis sets
sig_angio <- angio_gsea$pathway[angio_gsea$padj < 0.05]
if (length(sig_angio) > 0) {
  pdf(file.path(data_dir, "GSEA_angiogenesis_enrichment_plots.pdf"), width = 8, height = 5)
  for (pw in sig_angio) {
    p <- plotEnrichment(all_genesets[[pw]], ranked) +
      labs(title = gsub("HALLMARK_|KEGG_|REACTOME_", "", pw),
           subtitle = sprintf("NES = %.2f, padj = %.3f",
                              angio_gsea$NES[angio_gsea$pathway == pw],
                              angio_gsea$padj[angio_gsea$pathway == pw]))
    print(p)
  }
  dev.off()
}

# Core angiogenesis gene heatmap across conditions
core_angio_genes <- unique(unlist(all_genesets[angio_sets]))
core_angio_genes <- core_angio_genes[core_angio_genes %in% rownames(tumour)]

avg_expr <- AverageExpression(tumour_met_prim, features = core_angio_genes,
                               group.by = "sample", layer = "data")$RNA
# Z-score across samples
avg_z <- t(scale(t(avg_expr)))
# Keep top variable genes
top_angio <- names(sort(apply(avg_z, 1, var), decreasing = TRUE))[1:50]

col_annot <- HeatmapAnnotation(
  condition = pseudobulk$condition,
  col = list(condition = c("LPT_MET" = "#E63946", "Primary" = "#457B9D"))
)

pdf(file.path(data_dir, "heatmap_angiogenesis_genes.pdf"), width = 10, height = 14)
Heatmap(avg_z[top_angio, ],
        name            = "z-score",
        top_annotation  = col_annot,
        show_column_names = TRUE,
        row_names_gp    = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 9),
        cluster_rows    = TRUE,
        cluster_columns = TRUE,
        col             = circlize::colorRamp2(c(-2, 0, 2), c("#457B9D", "white", "#E63946")))
dev.off()

# Key angiogenesis genes — violin plots
key_angio_genes <- c("Vegfa", "Vegfb", "Vegfc", "Vegfd", "Kdr", "Flt1", "Flt4",
                     "Pdgfa", "Pdgfb", "Angpt1", "Angpt2", "Tie1", "Tek",
                     "Hif1a", "Epas1", "Dll4", "Notch1", "Notch4", "Nrp1", "Nrp2",
                     "Large1")
key_angio_genes <- key_angio_genes[key_angio_genes %in% rownames(tumour)]

pdf(file.path(data_dir, "violin_key_angiogenesis_LARGE1.pdf"), width = 16, height = 12)
VlnPlot(tumour_met_prim, features = key_angio_genes, group.by = "condition",
        cols   = c("LPT_MET" = "#E63946", "Primary" = "#457B9D"),
        pt.size = 0, ncol = 5) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ── 12. Module scoring: angiogenesis programme ────────────────────────────────
# Score each tumour cell for its angiogenesis programme activity
hallmark_angio <- all_genesets[["HALLMARK_ANGIOGENESIS"]]
hallmark_angio <- hallmark_angio[hallmark_angio %in% rownames(tumour)]

tumour <- AddModuleScore(tumour, features = list(hallmark_angio),
                          name = "angio_score", ctrl = 100)
names(tumour@meta.data)[grep("angio_score1", names(tumour@meta.data))] <- "angio_score"

pdf(file.path(data_dir, "angiogenesis_module_score.pdf"), width = 14, height = 5)
p1 <- FeaturePlot(tumour, features = "angio_score",
                  cols = c("grey90", "#E63946"), min.cutoff = "q5", max.cutoff = "q95") +
  ggtitle("Angiogenesis module score")
p2 <- VlnPlot(tumour_met_prim, features = "angio_score", group.by = "condition",
               cols = c("LPT_MET" = "#E63946", "Primary" = "#457B9D"), pt.size = 0) +
  ggtitle("Angiogenesis score: MET vs Primary") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0)
print(p1 | p2)
dev.off()

# ── 13. Save Seurat objects ───────────────────────────────────────────────────
# seurat_all_cells_integrated.rds was already saved in step 7, before `merged`
# was freed from memory.
saveRDS(tumour, file.path(data_dir, "seurat_tumour_cells.rds"))

message("\n=== Analysis complete ===")
message("Output files saved to: ", data_dir)
message("Key outputs:")
message("  - qc_violin.pdf")
message("  - umap_all_cells.pdf")
message("  - umap_tumour_classification.pdf")
message("  - umap_tumour_cells.pdf")
message("  - LARGE1_expression.pdf")
message("  - volcano_MET_vs_Primary.pdf")
message("  - GSEA_top_pathways.pdf")
message("  - GSEA_angiogenesis_enrichment_plots.pdf")
message("  - heatmap_angiogenesis_genes.pdf")
message("  - violin_key_angiogenesis_LARGE1.pdf")
message("  - angiogenesis_module_score.pdf")
message("  - DE_LPT_MET_vs_Primary.csv")
message("  - GSEA_MET_vs_Primary.csv")
message("  - GSEA_angiogenesis_pathways.csv")
message("  - LARGE1_per_patient.csv")
message("  - seurat_all_cells_integrated.rds")
message("  - seurat_tumour_cells.rds")
```

---

## Appendix B: Cluster labelling script (`R/downstream_analysis/label_clusters.R`)

Complete, unabridged, current version.

```r
# ── Manual cluster annotation (all-cells object) ─────────────────────────
# Labels derived from FindAllMarkers() output in tme_composition_analysis.R.
# Apply this after loading seurat_all_cells_integrated.rds.

library(Seurat)
library(dplyr)
library(ggplot2)

merged <- readRDS("seurat_all_cells_integrated.rds")

cluster_labels <- c(
  "0"  = "Tumour",
  "1"  = "Tumour",
  "2"  = "Smooth muscle",
  "3"  = "Tumour (S/G2M-phase)",
  "4"  = "Tumour (S/G2M-phase)",
  "5"  = "Pericytes",
  "6"  = "Tumour (S/G2M-phase)",
  "7"  = "Oligodendrocyte-like/glia",
  "8"  = "BBB endothelial",
  "9"  = "Schwann cells",
  "10" = "Fibroblasts",
  "11" = "Inflammatory macrophages/monocytes",
  "12" = "Microglia/macrophages",
  "13" = "Schwann cells",
  "14" = "Sex-linked/technical",
  "15" = "Tumour (stress/IEG state)",
  "16" = "Astrocytes",
  "17" = "Neutrophils",
  "18" = "T cells / NK cells",
  "19" = "Endothelial (Pecam1/Sox17-high)",
  "20" = "Homeostatic microglia",
  "21" = "Endothelial (Kdr/Mmrn2-high)",
  "22" = "Low-quality"
)

merged$cell_type <- factor(unname(cluster_labels[as.character(merged$seurat_clusters)]),
                            levels = unique(cluster_labels))

# Broad tumour vs TME flag, useful as a sanity check against CopyKAT's is_tumour
merged$compartment <- ifelse(grepl("^Tumour", merged$cell_type), "Tumour", "TME")

pdf("umap_all_cells_labelled.pdf", width = 10, height = 7)
print(DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE) +
        ggtitle("All cells: annotated clusters"))
print(DimPlot(merged, group.by = "compartment") +
        ggtitle("All cells: tumour vs TME"))
dev.off()

saveRDS(merged, "seurat_all_cells_integrated_labelled.rds")

# ── Composition by labelled cell type, MET vs Primary (paired by patient) ──
props <- merged@meta.data %>%
  filter(condition %in% c("LPT_MET", "Primary")) %>%
  count(patient, condition, cell_type) %>%
  group_by(patient, condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_wide <- props %>%
  select(patient, condition, cell_type, prop) %>%
  tidyr::pivot_wider(names_from = condition, values_from = prop, values_fill = 0)

cell_type_prop_test <- prop_wide %>%
  group_by(cell_type) %>%
  summarise(
    mean_MET     = mean(LPT_MET, na.rm = TRUE),
    mean_Primary = mean(Primary, na.rm = TRUE),
    p_value      = tryCatch(wilcox.test(LPT_MET, Primary, paired = TRUE)$p.value,
                             error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  arrange(p_value)

print(cell_type_prop_test, n = 30)
write.csv(cell_type_prop_test, "cell_type_proportion_test_MET_vs_Primary.csv", row.names = FALSE)

ggplot(props, aes(x = interaction(patient, condition), y = prop, fill = cell_type)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Patient / Condition", y = "Proportion of cells", fill = "Cell type",
       title = "Cell composition by patient and condition (labelled)")
ggsave("cell_type_composition_barplot.pdf", width = 11, height = 6)
```
