import os
import json
import pandas as pd
import numpy as np

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")

with open(f"{OUT_DIR}/paired_data.json") as f:
    paired = json.load(f)
stats_df = pd.read_csv(f"{OUT_DIR}/gene_stats_by_subgroup.csv")
gsea_df = pd.read_csv(f"{OUT_DIR}/gsea_angiogenesis_results.csv")

genes = ["LRG1", "CD74", "MIF", "EMILIN3", "ENG", "ITGB1", "PTK2"]
gene_label = {"PTK2": "PTK2 (FAK)", "ENG": "ENG (Endoglin)"}
group_order = ["group_3", "group_4", "shh"]
group_label = {"group_3": "Group 3 MB", "group_4": "Group 4 MB", "shh": "SHH-MB"}

PANEL_W, PANEL_H = 200, 190
PAD_L, PAD_R, PAD_T, PAD_B = 40, 14, 28, 34
INNER_W = PANEL_W - PAD_L - PAD_R
INNER_H = PANEL_H - PAD_T - PAD_B
X0, X1 = PAD_L + 20, PAD_L + INNER_W - 20


def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "n/a"
    if p < 0.001:
        return "p<0.001"
    return f"p={p:.3f}"


def make_panel(gene, grp, highlight=False):
    d = paired[grp][gene]
    prim = np.array(d["primary"])
    recur = np.array(d["relapse"])
    n = len(prim)
    row = stats_df[(stats_df["gene"] == gene) & (stats_df["group"] == group_label[grp])].iloc[0]
    pval = row["wilcoxon_p"]

    all_vals = np.concatenate([prim, recur])
    ymin, ymax = all_vals.min(), all_vals.max()
    span = ymax - ymin
    if span == 0:
        span = 1
    ymin -= span * 0.15
    ymax += span * 0.15

    def y(v):
        return PAD_T + INNER_H - (v - ymin) / (ymax - ymin) * INNER_H

    lines = []
    for pv, rv in zip(prim, recur):
        yp, yr = y(pv), y(rv)
        up = yr < yp
        color = "var(--series-2)" if up else "var(--series-1-muted)"
        lines.append(
            f'<line x1="{X0}" y1="{yp:.1f}" x2="{X1}" y2="{yr:.1f}" '
            f'stroke="{color}" stroke-width="1.4" opacity="0.75"/>'
        )
    dots = []
    for pv in prim:
        dots.append(f'<circle cx="{X0}" cy="{y(pv):.1f}" r="3.4" fill="var(--series-1)"/>')
    for rv in recur:
        dots.append(f'<circle cx="{X1}" cy="{y(rv):.1f}" r="3.4" fill="var(--series-2)"/>')

    # y-axis ticks (min/max)
    axis = (
        f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+INNER_H}" stroke="var(--axis-line)" stroke-width="1"/>'
        f'<text x="{PAD_L-6}" y="{PAD_T+4}" text-anchor="end" class="tick">{ymax:.1f}</text>'
        f'<text x="{PAD_L-6}" y="{PAD_T+INNER_H+2}" text-anchor="end" class="tick">{ymin:.1f}</text>'
    )
    xlabels = (
        f'<text x="{X0}" y="{PAD_T+INNER_H+16}" text-anchor="middle" class="tick">Primary</text>'
        f'<text x="{X1}" y="{PAD_T+INNER_H+16}" text-anchor="middle" class="tick">Relapse</text>'
    )
    sig = "sig" if (not np.isnan(pval) and pval < 0.05) else "nsig" if not np.isnan(pval) else ""
    title = f'<text x="{PANEL_W/2}" y="14" text-anchor="middle" class="panel-title">{group_label[grp]} (n={n})</text>'
    pnote = f'<text x="{PANEL_W/2}" y="{PANEL_H-4}" text-anchor="middle" class="pnote {sig}">{fmt_p(pval)}</text>'

    cls = "panel highlight" if highlight else "panel"
    return (
        f'<svg class="{cls}" width="{PANEL_W}" height="{PANEL_H}" viewBox="0 0 {PANEL_W} {PANEL_H}">'
        f'{title}{axis}{xlabels}{"".join(lines)}{"".join(dots)}{pnote}</svg>'
    )


rows_html = []
for gene in genes + ["angio_score"]:
    label = "Composite angiogenesis-gene score (mean z-score)" if gene == "angio_score" else f'{gene_label.get(gene, gene)}'
    panels = "".join(
        make_panel(gene, grp, highlight=(grp == "group_3"))
        for grp in group_order
    )
    caveat = ""
    if gene == "MIF":
        caveat = '<div class="caveat">Caveat: raw counts for MIF are near zero in almost all samples (uniquely-mapped-reads quantification likely undercounts MIF due to processed pseudogenes) — treat this gene\'s result as unreliable.</div>'
    rows_html.append(
        f'<section class="gene-row"><h2>{label}</h2>{caveat}<div class="panel-grid">{panels}</div></section>'
    )


# --- GSEA (preranked, paired t-stat relapse-vs-primary) bar chart ---
GSEA_W, GSEA_H = 620, 260
GBAR_PAD_L, GBAR_PAD_R, GBAR_PAD_T, GBAR_PAD_B = 190, 60, 10, 30
gbar_inner_w = GSEA_W - GBAR_PAD_L - GBAR_PAD_R
term_order = ["HALLMARK_ANGIOGENESIS", "CURATED_ANGIOGENESIS_SUPPLEMENTARY", "USER_ANGIOGENESIS_PANEL"]
term_display = {
    "HALLMARK_ANGIOGENESIS": "MSigDB Hallmark Angiogenesis (verified, 36 genes)",
    "CURATED_ANGIOGENESIS_SUPPLEMENTARY": "Supplementary curated angiogenesis panel (~50 genes, lower confidence)",
    "USER_ANGIOGENESIS_PANEL": "Your 7-gene panel (LRG1/CD74/MIF/EMILIN3/ENG/ITGB1/PTK2)",
}
grp_display_order = ["Group 3 MB", "Group 4 MB", "SHH-MB"]
series_colors = ["var(--series-1)", "var(--series-2)", "var(--series-3)"]

max_nes = gsea_df["NES"].abs().max() * 1.15
row_h = 26
n_rows = len(term_order) * len(grp_display_order)
GSEA_H = GBAR_PAD_T + n_rows * row_h + GBAR_PAD_B + 30

def x_scale(v):
    return GBAR_PAD_L + (v / max_nes) * gbar_inner_w

bars = []
y = GBAR_PAD_T + 20
zero_x = x_scale(0)
for grp in grp_display_order:
    bars.append(f'<text x="4" y="{y - 6}" class="gsea-group-label">{grp}</text>')
    for i, term in enumerate(term_order):
        row = gsea_df[(gsea_df["group"] == grp) & (gsea_df["Term"] == term)]
        if row.empty:
            y += row_h
            continue
        r = row.iloc[0]
        nes = r["NES"]
        fdr = r["FDR q-val"]
        bx = x_scale(min(nes, 0))
        bw = abs(x_scale(nes) - zero_x)
        sig = fdr < 0.05
        opacity = "1" if sig else "0.45"
        bars.append(
            f'<rect x="{bx:.1f}" y="{y-9:.1f}" width="{bw:.1f}" height="14" rx="3" '
            f'fill="{series_colors[i]}" opacity="{opacity}"/>'
        )
        label = f'NES {nes:.2f}{"  *FDR<0.05" if sig else ""}'
        bars.append(f'<text x="{x_scale(nes)+6:.1f}" y="{y+1:.1f}" class="gsea-val">{label}</text>')
        y += row_h
    y += 6

legend_items = "".join(
    f'<div class="legend-item"><span class="dot" style="background:{series_colors[i]}"></span>{term_display[t]}</div>'
    for i, t in enumerate(term_order)
)

gsea_svg = (
    f'<svg width="{GSEA_W}" height="{GSEA_H}" viewBox="0 0 {GSEA_W} {GSEA_H}">'
    f'<line x1="{zero_x:.1f}" y1="{GBAR_PAD_T}" x2="{zero_x:.1f}" y2="{GSEA_H-GBAR_PAD_B}" stroke="var(--axis-line)" stroke-width="1"/>'
    f'{"".join(bars)}'
    f'<text x="{zero_x:.1f}" y="{GSEA_H-8}" text-anchor="middle" class="tick">NES = 0</text>'
    f'</svg>'
)

gsea_section = f'''
<section class="gene-row">
<h2>GSEA (preranked): are angiogenesis gene sets shifted at relapse?</h2>
<p class="note">Genes ranked per subgroup by paired t-statistic (relapse&minus;primary, across n patients &mdash; Group 3/4 switch-at-relapse pairs bucketed by their primary-tumor group). Positive NES = the gene set skews toward genes <b>up</b> at relapse. Bars at full opacity pass FDR&lt;0.05.</p>
<div class="overflow-x">{gsea_svg}</div>
<div class="legend">{legend_items}</div>
</section>
'''

table_rows = []
for _, r in stats_df.iterrows():
    sig_class = "sig" if r["wilcoxon_p"] < 0.05 else ""
    glabel = "Angiogenesis score" if r["gene"] == "angio_score" else gene_label.get(r["gene"], r["gene"])
    table_rows.append(
        f'<tr class="{sig_class}"><td>{r["group"]}</td><td>{glabel}</td><td>{r["n_pairs"]}</td>'
        f'<td>{r["median_primary"]:.2f}</td><td>{r["median_relapse"]:.2f}</td>'
        f'<td>{r["log2FC_relapse_vs_primary"]:+.2f}</td><td>{r["n_increased_at_relapse"]}/{r["n_decreased_at_relapse"]}</td>'
        f'<td>{fmt_p(r["wilcoxon_p"])}</td></tr>'
    )

html = f"""<!doctype html>
<title>MB angiogenesis genes: primary vs relapse</title>
<style>
.viz-root {{
  color-scheme: light;
  --surface-1: #fcfcfb;
  --surface-2: #f3f2ef;
  --text-primary: #0b0b0b;
  --text-secondary: #52514e;
  --text-muted: #86847c;
  --axis-line: #d8d6cf;
  --series-1: #2a78d6;
  --series-1-muted: #a9c6ec;
  --series-2: #eb6834;
  --series-3: #1baf7a;
  --border: #e4e2db;
}}
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .viz-root {{
    color-scheme: dark;
    --surface-1: #1a1a19;
    --surface-2: #232321;
    --text-primary: #ffffff;
    --text-secondary: #c3c2b7;
    --text-muted: #8c8a80;
    --axis-line: #3a3934;
    --series-1: #3987e5;
    --series-1-muted: #3a5578;
    --series-2: #d95926;
    --series-3: #199e70;
    --border: #33322d;
  }}
}}
:root[data-theme="dark"] .viz-root {{
  color-scheme: dark;
  --surface-1: #1a1a19;
  --surface-2: #232321;
  --text-primary: #ffffff;
  --text-secondary: #c3c2b7;
  --text-muted: #8c8a80;
  --axis-line: #3a3934;
  --series-1: #3987e5;
  --series-1-muted: #3a5578;
  --series-2: #d95926;
  --series-3: #199e70;
  --border: #33322d;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; }}
.viz-root {{
  background: var(--surface-1);
  color: var(--text-primary);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  padding: 24px 20px 60px;
}}
.wrap {{ max-width: 900px; margin: 0 auto; }}
h1 {{ font-size: 20px; margin: 0 0 4px; }}
.subtitle {{ color: var(--text-secondary); font-size: 13px; margin: 0 0 6px; }}
.legend {{ display: flex; gap: 18px; align-items: center; font-size: 12px; color: var(--text-secondary); margin: 14px 0 24px; }}
.legend-item {{ display: flex; align-items: center; gap: 6px; }}
.dot {{ width: 9px; height: 9px; border-radius: 50%; display: inline-block; }}
.gene-row {{ margin-bottom: 30px; padding-bottom: 22px; border-bottom: 1px solid var(--border); }}
.gene-row:last-of-type {{ border-bottom: none; }}
h2 {{ font-size: 15px; margin: 0 0 10px; }}
.caveat {{ font-size: 12px; color: var(--text-muted); margin-bottom: 8px; max-width: 640px; }}
.panel-grid {{ display: flex; gap: 10px; flex-wrap: wrap; overflow-x: auto; }}
.panel {{ background: var(--surface-2); border-radius: 8px; border: 1px solid var(--border); }}
.panel.highlight {{ border: 1.5px solid var(--series-2); }}
.panel-title {{ font-size: 10px; fill: var(--text-secondary); }}
.tick {{ font-size: 9px; fill: var(--text-muted); }}
.pnote {{ font-size: 10px; fill: var(--text-muted); }}
.pnote.sig {{ fill: var(--series-2); font-weight: 600; }}
table {{ border-collapse: collapse; width: 100%; font-size: 12px; margin-top: 8px; }}
th, td {{ text-align: left; padding: 5px 8px; border-bottom: 1px solid var(--border); }}
th {{ color: var(--text-secondary); font-weight: 600; }}
tr.sig td {{ color: var(--series-2); font-weight: 600; }}
.overflow-x {{ overflow-x: auto; }}
.gsea-group-label {{ font-size: 11px; font-weight: 600; fill: var(--text-primary); }}
.gsea-val {{ font-size: 10px; fill: var(--text-secondary); }}
.note {{ font-size: 12px; color: var(--text-secondary); max-width: 700px; margin: 6px 0 20px; line-height: 1.5; }}
</style>
<div class="viz-root">
<div class="wrap">
<h1>Angiogenesis-related genes: primary vs relapse medulloblastoma</h1>
<p class="subtitle">Okonechnikov et al. 2023 cohort (n=43 primary-relapse pairs) &middot; log2 CPM from raw RNA-seq counts &middot; paired by patient &middot; Group 3 MB highlighted as primary focus</p>
<div class="legend">
  <div class="legend-item"><span class="dot" style="background:var(--series-1)"></span>Primary</div>
  <div class="legend-item"><span class="dot" style="background:var(--series-2)"></span>Relapse</div>
  <div class="legend-item">Line color = direction of change (orange = up at relapse, blue = down)</div>
</div>
{"".join(rows_html)}
{gsea_section}
<h2>Summary table</h2>
<p class="note">Wilcoxon signed-rank test (paired, primary vs relapse) per gene per subgroup. Rows in orange: p&lt;0.05. Group 3 MB has only 5 pairs so significance is hard to reach there even for consistent trends &mdash; check the n_increased/n_decreased column for consistency of direction.</p>
<div class="overflow-x">
<table>
<tr><th>Subgroup</th><th>Gene</th><th>n pairs</th><th>Median primary (log2CPM)</th><th>Median relapse (log2CPM)</th><th>&Delta; (relapse&minus;primary)</th><th>Up/Down (pairs)</th><th>Wilcoxon p</th></tr>
{"".join(table_rows)}
</table>
</div>
</div>
</div>
"""

with open(f"{OUT_DIR}/report.html", "w") as f:
    f.write(html)
print("written", f"{OUT_DIR}/report.html", len(html))
