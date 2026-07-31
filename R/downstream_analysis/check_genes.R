# ── Reusable gene-checking function ──────────────────────────────────────────
# Define once, then call check_genes(c("Gene1","Gene2",...)) as many times as
# you like with different gene lists — no need to re-load the object each time.

library(Seurat)
library(ggplot2)
library(dplyr)

# Adjust path if running locally vs on Myriad
tumour <- readRDS("seurat_tumour_cells.rds")
DefaultAssay(tumour) <- "SCT"

check_genes <- function(genes, split_by_cluster = FALSE,
                         out_file = paste0("check_genes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")) {
  found <- genes[genes %in% rownames(tumour)]
  missing <- setdiff(genes, found)
  if (length(missing) > 0) message("Not found, skipping: ", paste(missing, collapse = ", "))
  if (length(found) == 0) stop("None of the requested genes were found.")

  # Explicit named PDF device — avoids relying on R's default device, which
  # (when running via Rterm with no graphics window available) silently falls
  # back to writing "Rplots.pdf" and can end up corrupted if never closed.
  pdf(out_file, width = 12, height = 8)
  on.exit(dev.off(), add = TRUE)   # guarantees the device closes even on error

  # Violin plots by condition
  print(
    VlnPlot(tumour, features = found, group.by = "condition", pt.size = 0,
            ncol = min(4, length(found))) &
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  # UMAP feature plots
  print(
    FeaturePlot(tumour, features = found, cols = c("grey90", "#E63946"),
                min.cutoff = "q5", max.cutoff = "q95",
                ncol = min(4, length(found)))
  )

  # Optional: split by cluster too
  if (split_by_cluster) {
    print(
      DotPlot(tumour, features = found, group.by = "seurat_clusters") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    )
  }

  message("Plots written to: ", out_file)

  # Per-sample pseudobulk average (printed to console, not part of the PDF)
  avg <- AverageExpression(tumour, features = found,
                            group.by = c("sample", "condition"),
                            assays = "SCT", layer = "data")$SCT
  print(avg)

  invisible(avg)
}

# ── Usage ─────────────────────────────────────────────────────────────────
# check_genes(c("Vegfa", "Kdr", "Large1"))
# check_genes(c("Mki67", "Top2a"), split_by_cluster = TRUE)
# check_genes(c("Large1"), out_file = "large1_check.pdf")   # custom filename
