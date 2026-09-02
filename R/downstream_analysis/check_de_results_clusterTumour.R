# ── Pull a specific gene list out of the cluster-defined-tumour DESeq2 results
# Runs against DE_LPT_MET_vs_Primary_clusterTumour.csv (from
# rerun_gsea_cluster_based_tumour.R - tumour defined by all-cells clusters
# 0/1/3/4/6/15, NOT the original CopyKAT-based DE_LPT_MET_vs_Primary.csv,
# which carries ~14% TME contamination). Use this file when checking genes
# you want to trust, not the original CopyKAT-based one.

data_dir <- "."  # change to your local folder holding the downloaded files
de <- read.csv(file.path(data_dir, "DE_LPT_MET_vs_Primary_clusterTumour.csv"), row.names = 1)

genes_of_interest <- c("Lrg1", "Mif", "Cd74", "Emilin3", "Eng", "Itgb1", "Ptk2",
                        "Pglyrp1", "Lum", "Thbd", "Kcnj8", "Col3a1", "Pdgfa",
                        "Pdgfra", "App", "Trem1")
# Note: "Endoglin" is gene symbol Eng in mouse; "Mif/cd74" checked as two
# separate genes (Mif and its receptor Cd74) since DE is per-gene.

result <- de[de$gene %in% genes_of_interest, c("gene", "avg_log2FC", "p_val", "p_val_adj")]
result <- result[order(result$p_val_adj), ]

missing <- setdiff(genes_of_interest, result$gene)
if (length(missing) > 0) {
  message("Not present in DE results (likely filtered for low expression, or gene symbol mismatch): ",
          paste(missing, collapse = ", "))
}

print(result, row.names = FALSE)
write.csv(result, file.path(data_dir, "gene_list_DESeq2_results_clusterTumour.csv"), row.names = FALSE)
