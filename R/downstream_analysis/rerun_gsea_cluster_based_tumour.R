# ── Re-run pseudobulk DE + GSEA using a cluster-defined tumour population ────
# Replaces CopyKAT's `is_tumour` call entirely with a manual cluster-based
# tumour definition (clusters 0, 1, 3, 4, 6, 15 from the all-cells UMAP,
# resolution 0.5), motivated by CopyKAT's ~14% TME contamination rate found
# via marker-gene cross-check earlier in this project.
#
# Cluster 14 was deliberately excluded: its real top-10 markers (Tsix,
# Gm26917, Gria2, Miat, Gm42418, Xist, Kcnq1ot1, Sacs, C130071C03Rik, Nav2)
# contain no tumour-lineage marker (no Neurod1/Atoh1/Barhl1/Cntn2), are
# dominated by sex-chromosome/imprinted lncRNAs, and include Gm42418 - a
# high-expression but non-specific rRNA-repeat transcript that's a common
# technical-artifact signature in mouse scRNA-seq, not a real marker.
#
# Runs entirely locally - no HPC/Myriad needed. Requires: Seurat, dplyr,
# ggplot2, DESeq2, fgsea, msigdbr installed locally.

library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(fgsea)
library(msigdbr)

# ── 0. Load the all-cells object (adjust path to wherever you downloaded it) ─
data_dir <- "."  # change to your local folder holding the downloaded files
merged <- readRDS(file.path(data_dir, "seurat_all_cells_integrated.rds"))

# ── 1. Define tumour cells by cluster, not by CopyKAT ────────────────────────
tumour_clusters <- c(0, 1, 3, 4, 6, 15)
tumour <- subset(merged, subset = seurat_clusters %in% tumour_clusters)
rm(merged); gc()

message("Cluster-defined tumour object: ", ncol(tumour), " cells across ",
        length(unique(tumour$sample)), " samples")
print(table(tumour$seurat_clusters, tumour$condition))

# ── 2. Pseudobulk DE: LPT_MET vs Primary ─────────────────────────────────────
DefaultAssay(tumour) <- "RNA"
tumour_met_prim <- subset(tumour, subset = condition %in% c("LPT_MET", "Primary"))

pseudobulk <- AggregateExpression(tumour_met_prim,
  group.by = c("sample", "condition", "patient"),
  assays    = "RNA",
  return.seurat = TRUE
)

Idents(pseudobulk) <- "condition"
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
de_results$gene <- rownames(de_results)
de_results <- de_results[order(de_results$avg_log2FC, decreasing = TRUE), ]
write.csv(de_results, file.path(data_dir, "DE_LPT_MET_vs_Primary_clusterTumour.csv"), row.names = TRUE)

# ── 3. GSEA (same defensive gene-set loading as the main pipeline) ──────────
de_results$rank_metric <- de_results$avg_log2FC * -log10(de_results$p_val + 1e-300)
ranked <- setNames(de_results$rank_metric, de_results$gene)
ranked <- sort(ranked, decreasing = TRUE)

load_geneset <- function(...) {
  tryCatch({
    msigdbr(...) %>%
      dplyr::select(gs_name, gene_symbol) %>%
      split(x = .$gene_symbol, f = .$gs_name)
  }, error = function(e) {
    message("Skipping gene set (", paste(c(...), collapse = "/"), "): ", conditionMessage(e))
    list()
  })
}

geneset_groups <- list(
  hallmark = list(
    list(species = "Mus musculus", db_species = "MM", collection = "MH"),
    list(species = "Mus musculus", collection = "H")
  ),
  kegg = list(
    list(species = "Mus musculus", db_species = "MM", collection = "M2", subcollection = "CP:KEGG_LEGACY"),
    list(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  ),
  reactome = list(
    list(species = "Mus musculus", db_species = "MM", collection = "M2", subcollection = "CP:REACTOME"),
    list(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
  )
)

all_genesets <- list()
for (group_name in names(geneset_groups)) {
  for (spec in geneset_groups[[group_name]]) {
    gs <- do.call(load_geneset, spec)
    if (length(gs) > 0) {
      all_genesets <- modifyList(all_genesets, gs)
      break
    }
  }
}
if (length(all_genesets) == 0) stop("No MSigDB gene sets could be loaded.")
message(length(all_genesets), " MSigDB gene sets loaded for GSEA")

set.seed(42)
gsea_res <- fgsea(pathways = all_genesets, stats = ranked,
                   minSize = 10, maxSize = 500, nPermSimple = 10000)
gsea_res <- gsea_res[order(gsea_res$NES, decreasing = TRUE), ]
write.csv(gsea_res[, -8], file.path(data_dir, "GSEA_MET_vs_Primary_clusterTumour.csv"), row.names = FALSE)

# ── 4. Angiogenesis focus (same keyword filter as the main pipeline) ────────
angio_sets <- grep("ANGIOGEN|VEGF|NOTCH|HIF|HYPOXIA|VESSEL|VASCULO",
                   names(all_genesets), value = TRUE, ignore.case = TRUE)
angio_gsea <- gsea_res[gsea_res$pathway %in% angio_sets, ]
write.csv(angio_gsea[, -8], file.path(data_dir, "GSEA_angiogenesis_pathways_clusterTumour.csv"), row.names = FALSE)

print(as.data.frame(angio_gsea[, c("pathway", "NES", "padj")]))

# ── 5. Compare against the CopyKAT-based results, if you have that CSV too ──
copykat_de_path <- file.path(data_dir, "DE_LPT_MET_vs_Primary.csv")
if (file.exists(copykat_de_path)) {
  copykat_de <- read.csv(copykat_de_path)
  merged_de <- merge(de_results, copykat_de, by = "gene", suffixes = c("_clusterTumour", "_copykat"))
  write.csv(merged_de, file.path(data_dir, "DE_comparison_clusterTumour_vs_copykat.csv"), row.names = FALSE)
  message("Wrote side-by-side comparison to DE_comparison_clusterTumour_vs_copykat.csv")
}

# ── 6. Genome-wide top-pathways plot (not restricted to angiogenesis) ───────
# Mirrors GSEA_top_pathways.pdf from the main Myriad pipeline, but for this
# cluster-defined-tumour result. Useful for presenting the whole-transcriptome
# signal (e.g. HALLMARK_TNFA_SIGNALING_VIA_NFKB, the translation/ribosome/
# OXPHOS cluster) rather than only the angiogenesis-filtered subset.
top_paths <- rbind(
  head(gsea_res[gsea_res$padj < 0.05 & gsea_res$NES > 0, ], 15),
  head(gsea_res[gsea_res$padj < 0.05 & gsea_res$NES < 0, ], 15)
)
if (nrow(top_paths) > 0) {
  top_paths$pathway_short <- gsub("HALLMARK_|KEGG_|REACTOME_", "", top_paths$pathway)
  top_paths$pathway_short <- substr(top_paths$pathway_short, 1, 55)
  top_paths$direction <- ifelse(top_paths$NES > 0, "Up in MET", "Up in Primary")

  pdf(file.path(data_dir, "GSEA_top_pathways_clusterTumour.pdf"), width = 11, height = 10)
  print(
    ggplot(top_paths, aes(x = NES, y = reorder(pathway_short, NES), fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c("Up in MET" = "#E63946", "Up in Primary" = "#457B9D")) +
      geom_vline(xintercept = 0, colour = "black") +
      labs(x = "Normalised Enrichment Score", y = NULL,
           title = "GSEA: LPT_MET vs Primary (cluster-defined tumour, all significant pathways)") +
      theme_classic() + theme(legend.title = element_blank())
  )
  dev.off()
  message("Wrote GSEA_top_pathways_clusterTumour.pdf (", nrow(top_paths), " significant pathways, padj < 0.05)")
} else {
  message("No pathways reached padj < 0.05 - skipping top-pathways plot")
}

# ── 7. Sequencing-depth QC: rule out a technical confound before trusting ───
# the translation/ribosome/OXPHOS GSEA signal. A large, coordinated shift in
# ribosomal/mitochondrial-translation gene sets is a classic signature of a
# per-sample sequencing-depth or RNA-quality difference, not necessarily real
# biology - check nCount_RNA/nFeature_RNA by condition before presenting it.
depth_by_sample <- tumour_met_prim@meta.data %>%
  group_by(sample, condition, patient) %>%
  summarise(
    median_nCount   = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    n_cells         = n(),
    .groups = "drop"
  )
print(as.data.frame(depth_by_sample))
write.csv(depth_by_sample, file.path(data_dir, "sequencing_depth_by_sample_clusterTumour.csv"), row.names = FALSE)

depth_test <- function(metric) {
  wide <- depth_by_sample %>%
    select(patient, condition, !!metric) %>%
    tidyr::pivot_wider(names_from = condition, values_from = !!metric)
  tryCatch(wilcox.test(wide$LPT_MET, wide$Primary, paired = TRUE)$p.value,
           error = function(e) NA_real_)
}
message("Paired Wilcoxon, MET vs Primary - median nCount_RNA per sample: p = ",
        round(depth_test("median_nCount"), 4))
message("Paired Wilcoxon, MET vs Primary - median nFeature_RNA per sample: p = ",
        round(depth_test("median_nFeature"), 4))

pdf(file.path(data_dir, "sequencing_depth_by_condition_clusterTumour.pdf"), width = 8, height = 5)
print(
  ggplot(depth_by_sample, aes(x = condition, y = median_nCount, fill = condition)) +
    geom_boxplot() + geom_point(position = position_jitter(width = 0.1)) +
    theme_classic() + labs(title = "Median nCount_RNA per sample, by condition")
)
print(
  ggplot(depth_by_sample, aes(x = condition, y = median_nFeature, fill = condition)) +
    geom_boxplot() + geom_point(position = position_jitter(width = 0.1)) +
    theme_classic() + labs(title = "Median nFeature_RNA per sample, by condition")
)
dev.off()
message("Wrote sequencing_depth_by_condition_clusterTumour.pdf - check whether MET/Primary differ",
        " in depth before trusting the translation/ribosome/OXPHOS GSEA signal")
