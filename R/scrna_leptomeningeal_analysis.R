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
# msigdbr >= v10 renamed category/subcategory to collection/subcollection,
# and renamed the KEGG legacy subcollection code from "CP:KEGG" to
# "CP:KEGG_LEGACY" - the old names now throw "Unknown subcollection".
# With db_species = "MM" (mouse-native gene sets, rather than human sets
# ortholog-mapped to mouse), collection codes also change prefix: Hallmark
# is "MH" instead of "H", and the C2 curated collection is "M2" instead of
# "C2" - using the human-prefixed codes here throws "Unknown collection".
hallmark <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

kegg <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "M2", subcollection = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

reactome <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "M2", subcollection = "CP:REACTOME") %>%
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
