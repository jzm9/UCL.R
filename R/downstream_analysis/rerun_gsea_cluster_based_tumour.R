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

print(angio_gsea[, c("pathway", "NES", "padj")], n = 30)

# ── 5. Compare against the CopyKAT-based results, if you have that CSV too ──
copykat_de_path <- file.path(data_dir, "DE_LPT_MET_vs_Primary.csv")
if (file.exists(copykat_de_path)) {
  copykat_de <- read.csv(copykat_de_path)
  merged_de <- merge(de_results, copykat_de, by = "gene", suffixes = c("_clusterTumour", "_copykat"))
  write.csv(merged_de, file.path(data_dir, "DE_comparison_clusterTumour_vs_copykat.csv"), row.names = FALSE)
  message("Wrote side-by-side comparison to DE_comparison_clusterTumour_vs_copykat.csv")
}
