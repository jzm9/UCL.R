# Secondary vs primary tumour — GO:BP enrichment

Extends the earlier angiogenesis/Lrg1-only check to the full GO:BP (Gene
Ontology Biological Process) collection from MSigDB, so we can see which
processes actually top the ranked gene list in either direction, not just
angiogenesis.

- `data/secondary_vs_primary.csv` — edgeR-style DE table (contrast:
  secondary vs primary; `logFC > 0` = higher in secondary).
- `scripts/gobp_gsea.R` — loads the DE table, checks Lrg1, runs `fgsea`
  against every `GO:BP` gene set from `msigdbr` (mouse), and reports the
  top enriched terms in each direction plus where angiogenesis/vascular
  terms specifically rank.

Run from `scripts/` (paths in the script are relative to that directory):

```r
setwd("analyses/secondary-vs-primary/scripts")
source("gobp_gsea.R")
```

Requires `data.table`, `fgsea`, `msigdbr` (same packages used in the
original angiogenesis script). Note: `msigdbr` >= 10 needs the `msigdbdf`
backend package installed, and its API changed from
`category`/`subcategory` to `collection`/`subcollection` — the script
tries the new form first and falls back to the old one.

Outputs `results/gobp_gsea_full.csv` (every tested GO:BP term, NES/p/padj,
sorted by padj).
