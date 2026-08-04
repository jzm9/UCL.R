# =============================================================================
# scRNA-seq analysis: mouse leptomeningeal metastases vs primary tumours vs
# non-tumour leptomeninges (GSE264326, 11 samples)
#
# Groups:
#   LPT_MET_1-4  : leptomeningeal metastases   (GSM8216129-32ish)
#   Primary_1-4  : primary tumours (matched patients)
#   LPT_WT_1-3   : non-tumour leptomeninges (unpaired baseline)
#
# Run on UCL Myriad HPC via SGE (see run_scrna_myriad.sh). Designed to run
# end-to-end in a single batch job: QC -> per-sample SCTransform -> merge ->
# Harmony integration -> clustering -> CopyKAT tumour calling -> tumour
# subset re-analysis -> LARGE1 test -> pseudobulk DESeq2 -> GSEA -> targeted
# angiogenesis follow-up.
#
# NOTE ON THE ALL-CELLS OBJECT: seurat_all_cells_integrated.rds is saved
# straight after integration/clustering, BEFORE CopyKAT tumour filtering.
# CopyKAT-based tumour/non-tumour calling was later found to leave ~14% TME
# (non-tumour) cells misclassified as "tumour" in the filtered object, so any
# analysis of tumour microenvironment (TME) composition should use the
# all-cells object with is_tumour as metadata, not the filtered subset.
# =============================================================================

## ---- 0. Packages ------------------------------------------------------
required_pkgs <- c("Seurat", "SeuratObject", "harmony", "dplyr", "ggplot2",
                    "Matrix", "DESeq2", "fgsea", "msigdbr", "copykat")
# (clusterProfiler intentionally omitted - not needed, caused install issues)
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Missing package: ", p, " - install before running.")
  }
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

data_dir <- "/home/sejkmor/Scratch/taylordataset/taylor25"
setwd(data_dir)

## ---- 1. Sample metadata ------------------------------------------------
# 11 samples total: 4 metastasis / 4 primary (patient-matched) + 3 non-tumour
# leptomeninges baseline (unpaired).
sample_table <- data.frame(
  sample    = c("LPT_MET_1", "LPT_MET_2", "LPT_MET_3", "LPT_MET_4",
                "Primary_1", "Primary_2", "Primary_3", "Primary_4",
                "LPT_WT_1", "LPT_WT_2", "LPT_WT_3"),
  condition = c(rep("LPT_MET", 4), rep("Primary", 4), rep("LPT_WT", 3)),
  patient   = c(1, 2, 3, 4, 1, 2, 3, 4, NA, NA, NA),
  gsm       = c("GSM8216129", "GSM8216130", "GSM8216131", "GSM8216132",
                "GSM8216133", "GSM8216134", "GSM8216135", "GSM8216136",
                "GSM8216137", "GSM8216138", "GSM8216139"),
  stringsAsFactors = FALSE
)

## ---- 2. Load + QC filter each sample ------------------------------------
qc_summary <- list()
seurat_list <- list()

for (i in seq_len(nrow(sample_table))) {
  s <- sample_table$sample[i]
  counts <- Read10X(data.dir = file.path(data_dir, "raw", s))
  obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 3, min.features = 200)
  obj$sample    <- s
  obj$condition <- sample_table$condition[i]
  obj$patient   <- sample_table$patient[i]

  obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")  # mouse mito prefix

  qc_summary[[s]] <- data.frame(
    sample = s, n_cells_pre = ncol(obj),
    median_nFeature = median(obj$nFeature_RNA),
    median_nCount   = median(obj$nCount_RNA),
    median_pct_mt   = median(obj$pct_mt)
  )

  # QC thresholds
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                              nCount_RNA > 500 & pct_mt < 20)

  seurat_list[[s]] <- obj
}

write.csv(do.call(rbind, qc_summary), file.path(data_dir, "qc_summary_pre_filter.csv"), row.names = FALSE)

pdf(file.path(data_dir, "qc_violin.pdf"), width = 14, height = 6)
print(VlnPlot(merge(seurat_list[[1]], seurat_list[-1]),
              features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
              group.by = "sample", pt.size = 0, ncol = 3))
dev.off()

## ---- 3. Per-sample SCTransform normalisation ----------------------------
# Cap at 5000 cells/sample for tractable memory use on the cluster.
for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  if (ncol(obj) > 5000) obj <- subset(obj, cells = sample(colnames(obj), 5000))
  obj <- SCTransform(obj, vst.flavor = "v2", variable.features.n = 3000, verbose = FALSE)
  seurat_list[[s]] <- obj
}

## ---- 4. Merge + Harmony batch integration --------------------------------
merged <- merge(seurat_list[[1]], seurat_list[-1], merge.data = TRUE)
VariableFeatures(merged) <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- harmony::RunHarmony(merged, group.by.vars = "sample", max_iter = 20)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)

pdf(file.path(data_dir, "umap_all_cells.pdf"), width = 8, height = 6)
print(DimPlot(merged, group.by = "seurat_clusters", label = TRUE))
print(DimPlot(merged, group.by = "condition"))
dev.off()

# Save the full, unfiltered object BEFORE tumour subsetting. This is the
# object to use for any TME-composition analysis, since CopyKAT filtering
# below is imperfect (~14% TME contamination in the "tumour" call).
saveRDS(merged, file.path(data_dir, "seurat_all_cells_integrated.rds"))

## ---- 5. CopyKAT CNV-based tumour/non-tumour classification --------------
# genome = "mm10" is required for mouse data (default assumes human hg20).
library(copykat)
raw_counts <- GetAssayData(merged, assay = "RNA", layer = "counts")
copykat_res <- copykat(rawmat = as.matrix(raw_counts), id.type = "S",
                        genome = "mm10", sam.name = "leptomeningeal",
                        n.cores = as.integer(Sys.getenv("NSLOTS", "1")))

pred <- copykat_res$prediction
merged$is_tumour <- pred$copykat.pred[match(colnames(merged), pred$cell.names)]

pdf(file.path(data_dir, "umap_tumour_classification.pdf"), width = 8, height = 6)
print(DimPlot(merged, group.by = "is_tumour"))
dev.off()

## ---- 6. Tumour-only subset: re-normalise, re-integrate, re-cluster ------
tumour <- subset(merged, subset = is_tumour == "aneuploid")  # CopyKAT "tumour" call
DefaultAssay(tumour) <- "RNA"
tumour <- SCTransform(tumour, vst.flavor = "v2", variable.features.n = 3000, verbose = FALSE)
tumour <- RunPCA(tumour, npcs = 30, verbose = FALSE)
tumour <- harmony::RunHarmony(tumour, group.by.vars = "sample", max_iter = 20)
tumour <- RunUMAP(tumour, reduction = "harmony", dims = 1:30)
tumour <- FindNeighbors(tumour, reduction = "harmony", dims = 1:30)
tumour <- FindClusters(tumour, resolution = 0.4)

pdf(file.path(data_dir, "umap_tumour_cells.pdf"), width = 8, height = 6)
print(DimPlot(tumour, group.by = "seurat_clusters", label = TRUE))
print(DimPlot(tumour, group.by = "condition"))
dev.off()

## ---- 7. LARGE1 expression: LPT_MET vs Primary, paired by patient --------
DefaultAssay(tumour) <- "SCT"
tumour_met_prim <- subset(tumour, subset = condition %in% c("LPT_MET", "Primary"))

pdf(file.path(data_dir, "LARGE1_expression.pdf"), width = 6, height = 5)
print(VlnPlot(tumour_met_prim, features = "Large1", group.by = "condition", pt.size = 0))
dev.off()

large1_per_patient <- AverageExpression(tumour_met_prim, features = "Large1",
                                         group.by = c("patient", "condition"),
                                         assays = "SCT", layer = "data")$SCT
write.csv(large1_per_patient, file.path(data_dir, "LARGE1_per_patient.csv"))

large1_df <- FetchData(tumour_met_prim, vars = c("Large1", "patient", "condition")) %>%
  group_by(patient, condition) %>% summarise(mean_expr = mean(Large1), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_expr)
large1_wilcox <- wilcox.test(large1_df$LPT_MET, large1_df$Primary, paired = TRUE)
print(large1_wilcox)

## ---- 8. Pseudobulk DESeq2: LPT_MET vs Primary (tumour cells) ------------
tumour_met_prim <- PrepSCTFindMarkers(tumour_met_prim)  # required: multiple per-sample SCT models

pb <- AggregateExpression(tumour_met_prim, assays = "RNA", group.by = c("sample", "condition"),
                           return.seurat = TRUE)

ident_met  <- grep("LPT_MET", colnames(pb), value = TRUE)
ident_prim <- grep("Primary", colnames(pb), value = TRUE)

de_results <- FindMarkers(pb, ident.1 = ident_met, ident.2 = ident_prim,
                           test.use = "DESeq2", slot = "counts")
de_results$gene <- rownames(de_results)
write.csv(de_results, file.path(data_dir, "DE_LPT_MET_vs_Primary.csv"), row.names = FALSE)

# Volcano plot
de_results$sig <- with(de_results, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
pdf(file.path(data_dir, "volcano_DE_MET_vs_Primary.pdf"), width = 7, height = 6)
print(
  ggplot(de_results, aes(avg_log2FC, -log10(p_val_adj), color = sig)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    theme_classic() + labs(title = "LPT_MET vs Primary (pseudobulk DESeq2)")
)
dev.off()

## ---- 9. GSEA on the pseudobulk DE ranks (angiogenesis focus) ------------
library(fgsea)
library(msigdbr)

# msigdbr >= v10 API uses collection/subcollection, not category/subcategory.
# The KEGG legacy gene set collection was renamed CP:KEGG_LEGACY.
msig <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
msig_hallmark <- msigdbr(species = "Mus musculus", collection = "H")
msig_all <- rbind(msig, msig_hallmark)
pathway_list <- split(msig_all$gene_symbol, msig_all$gs_name)

ranks <- de_results$avg_log2FC
names(ranks) <- de_results$gene
ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)

gsea_res <- fgsea(pathways = pathway_list, stats = ranks, minSize = 10, maxSize = 500)
gsea_angio <- gsea_res[grepl("ANGIOGENESIS|VEGF|VASCULAR", gsea_res$pathway, ignore.case = TRUE), ]
gsea_angio <- gsea_angio[order(gsea_angio$padj), ]
write.csv(gsea_angio, file.path(data_dir, "GSEA_angiogenesis_pathways.csv"), row.names = FALSE)

pdf(file.path(data_dir, "GSEA_angiogenesis_enrichment_plots.pdf"), width = 7, height = 5)
for (pw in gsea_angio$pathway[seq_len(min(5, nrow(gsea_angio)))]) {
  print(plotEnrichment(pathway_list[[pw]], ranks) + ggtitle(pw))
}
dev.off()

## ---- 10. Angiogenesis gene deep-dive (heatmap + violins) ------------------
angio_genes <- unique(unlist(pathway_list[grepl("ANGIOGENESIS", names(pathway_list))]))
angio_genes_present <- intersect(angio_genes, rownames(tumour[["SCT"]]))

avg_by_sample <- AverageExpression(tumour_met_prim, features = angio_genes_present,
                                    group.by = c("sample", "condition"),
                                    assays = "SCT", layer = "data")$SCT

# Drop genes with zero variance across samples (undetected, or identical in
# every sample) - scale="row" divides by SD, so these become NaN/Inf and
# make pheatmap render a blank page instead of erroring.
heatmap_mat <- log1p(as.matrix(avg_by_sample))
row_sd <- apply(heatmap_mat, 1, sd)
heatmap_mat <- heatmap_mat[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
message(nrow(heatmap_mat), " / ", nrow(avg_by_sample),
        " angiogenesis genes have non-zero variance and will be plotted")

pdf(file.path(data_dir, "heatmap_angiogenesis_genes.pdf"), width = 10, height = 12)
if (nrow(heatmap_mat) >= 2) {
  pheatmap::pheatmap(heatmap_mat, scale = "row", show_rownames = TRUE)
} else {
  message("Fewer than 2 genes with variance - skipping heatmap")
}
dev.off()

# Key genes, plotted individually with tryCatch since not all are present /
# not all have non-degenerate expression in every group (VlnPlot crashes on
# all-identical values otherwise).
key_angio_genes <- c("Vegfa", "Kdr", "Flt1", "Angpt2", "Tek", "Pecam1", "Large1")
pdf(file.path(data_dir, "violin_key_angiogenesis_LARGE1.pdf"), width = 12, height = 8)
for (g in key_angio_genes) {
  tryCatch({
    print(VlnPlot(tumour_met_prim, features = g, group.by = "condition", pt.size = 0))
  }, error = function(e) message("Skipped ", g, ": ", conditionMessage(e)))
}
dev.off()

## ---- 11. Manual angiogenesis module score ---------------------------------
# Written by hand instead of AddModuleScore() due to a Seurat/S4 ABI bug
# ([[<- defined for objects of type "S4" only for subclasses of environment)
# encountered on this Seurat/R build.
angio_mat <- GetAssayData(tumour, assay = "SCT", layer = "data")[angio_genes_present, ]
angio_mat_scaled <- t(scale(t(as.matrix(angio_mat))))
angio_mat_scaled[is.na(angio_mat_scaled)] <- 0
tumour$angio_score <- colMeans(angio_mat_scaled)
tumour_met_prim$angio_score <- tumour$angio_score[colnames(tumour_met_prim)]

pdf(file.path(data_dir, "angiogenesis_module_score.pdf"), width = 6, height = 5)
print(VlnPlot(tumour_met_prim, features = "angio_score", group.by = "condition", pt.size = 0))
dev.off()

## ---- 12. Save final tumour object -----------------------------------------
saveRDS(tumour, file.path(data_dir, "seurat_tumour_cells.rds"))

cat("Pipeline complete.\n")
