# MB primary vs. relapse: angiogenesis genes and pathways

Standalone analysis, unrelated to the `velocyto.R` package itself — kept in its
own folder here at the user's request rather than as a separate repo.

Source data: Okonechnikov et al. 2023, *Acta Neuropathologica Communications*
11:7 ("Comparison of transcriptome profiles between medulloblastoma primary
and recurrent tumors uncovers novel variance effects in relapses"). Raw RNA-seq
counts (55,765 genes x 86 samples, 43 primary-relapse pairs: 24 SHH-MB, 14
Group 4 MB, 5 Group 3 MB, including 2 pairs that switch Group 4 -> Group 3 at
relapse) plus sample metadata (status/subgroup/patient pairing/clinical
variables) are hosted on the R2 Genomics platform
("Tumor Medulloblastoma—Korshunov—86—rpkm—mbffpe", http://r2.amc.nl) and are
**not committed to this repo** (patient-level sequencing data from a
consented, ethics-approved cohort). `results/sample_metadata.csv` (parsed
clinical/molecular annotation only, no sequencing data) is included since it's
needed to reproduce the derived results below.

## Pipeline

1. `scripts/01_build_metadata.py` — parses the R2 datagrabber metadata file
   into `results/sample_metadata.csv` (sample -> status/subgroup/tumorid/age/
   PFS/etc.), matched against the raw counts file's sample IDs.
2. `scripts/02_load_counts.py` — loads raw counts, converts to log2(CPM+1)
   (library-size normalized; simpler than the paper's RPKM but adequate for
   paired same-gene comparisons).
3. `scripts/03_analyze_genes.py` — paired Wilcoxon signed-rank test (primary
   vs. relapse, same patient) per gene per subgroup, for a set of user-chosen
   genes (LRG1, CD74, MIF, EMILIN3, ENG, ITGB1, PTK2) plus a composite
   z-score across them -> `results/gene_stats_by_subgroup.csv`.
4. `scripts/05_gsea.py` — preranked GSEA (genes ranked per subgroup by paired
   t-statistic, relapse-vs-primary) against:
   - `HALLMARK_ANGIOGENESIS` (MSigDB Hallmark v2023.1, 36 genes — verified via
     matching SHA-450d222 across independent public mirrors of the official
     Broad Institute release, since this environment cannot reach
     gsea-msigdb.org or Enrichr directly)
   - a supplementary curated angiogenesis panel (~50 genes, lower confidence —
     hand-assembled, not verified against an official GO Biological Process
     term due to file-size limits on retrieval)
   - the user's 7-gene panel
   -> `results/gsea_angiogenesis_results.csv`
5. `scripts/04_make_report.py` — builds `results/report.html`, a self-contained
   interactive report with paired slope plots (primary -> relapse per patient)
   for each gene x subgroup, plus the GSEA bar chart. Also published as a
   Claude artifact.

## Key findings

- **ITGB1** is up at relapse consistently across all three subgroups
  (significant in SHH-MB, p<0.001; directionally consistent 5/5 in Group 3 MB
  and 8/12 in Group 4 MB).
- The **Hallmark Angiogenesis** gene set as a whole is significantly enriched
  toward relapse in **Group 4 MB** (NES 2.08, FDR<0.001) and **SHH-MB**
  (NES 1.53, FDR=0.038), with the same positive direction (not significant,
  small n=5) in **Group 3 MB**.
- **MIF** raw counts are near-zero in almost all samples — likely a
  quantification artifact (uniquely-mapped-reads-only counting undercounts MIF
  due to its processed pseudogenes) — treat that gene's result as unreliable.

## Reproducing

Point `scripts/01_build_metadata.py` and `scripts/02_load_counts.py` at your
local copies of the R2 datagrabber metadata file and the raw counts file, run
scripts 01 -> 05 in order (Python 3, `pandas numpy scipy gseapy`).
