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
2. Top-10 markers per cluster (by `avg_log2FC`) inspected manually.
3. Identities assigned by matching marker genes to literature/canonical
   cell-type markers, by eye, cluster by cluster (`R/downstream_analysis/label_clusters.R`).
4. Adjacent/overlapping identities (e.g. clusters 19 and 21, both initially
   labelled "Endothelial (vascular)") were revisited once their actual
   top-10 marker sets were checked: 19 is defined by `Pecam1`/`Sox17`
   (canonical pan-endothelial), 21 by `Kdr`/`Mmrn2`/`Bvht`/`Ch25h` — non-
   overlapping gene sets, split into separate labels rather than merged,
   since collapsing them would have hidden a real subpopulation distinction.

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
