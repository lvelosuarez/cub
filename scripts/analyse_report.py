#!/usr/bin/env python3
"""
analyse_report.py — TPOS QC HTML report generator.

Single public function: write_tpos_html().
No external dependencies beyond the Python stdlib and polars.
"""
from pathlib import Path
import datetime
import polars as pl


def write_tpos_html(
    tpos_run_sample_df: pl.DataFrame,
    comparison_df:      pl.DataFrame,
    qc_row:             pl.DataFrame,
    run_id:             int,
    out_path:           Path,
) -> None:
    """
    Write a self-contained HTML TPOS QC report.

    Parameters
    ----------
    tpos_run_sample_df : raw run_sample() output for the TPOS sample
    comparison_df      : expected vs observed organism table from run_tpos_qc()
    qc_row             : single-row DataFrame from run_tpos_qc()
    run_id             : integer run identifier
    out_path           : destination path for the HTML file
    """
    r = qc_row.row(0, named=True)
    qc_flag   = r["qc_flag"]
    rho       = r["spearman_rho"]
    lod       = r["lod_fraction"]
    must_rate = r["must_detect_rate"]
    n_det     = r["n_detected"]
    n_exp     = r["n_expected"]
    today     = datetime.date.today().isoformat()

    flag_color = {"PASS": "#22c55e", "WARN": "#f59e0b", "FAIL": "#ef4444"}.get(qc_flag, "#6b7280")

    # Unexpected organisms: non-TPOS taxa in the raw TPOS run_sample output
    TPOS_TAXIDS = frozenset([1639, 287, 1423, 4932, 562, 28901, 1613, 1351, 5207, 1280])
    top_unexpected = (
        tpos_run_sample_df
        .filter(~pl.col("group_taxid").is_in(list(TPOS_TAXIDS)))
        .sort("n_reads_total", descending=True)
        .head(10)
    )

    def _td(val, cls=""):
        return f'<td class="{cls}">{val}</td>'

    def _pct(v):
        return f"{v:.3%}" if v is not None else "—"

    def _fc(v):
        return f"{v:+.1f}" if v is not None else "n/d"

    # Build comparison rows
    cmp_rows_html = []
    for row in comparison_df.iter_rows(named=True):
        det    = row.get("observed_fraction", 0) or 0
        sym    = "&#10003;" if det > 0 else "&#10007;"
        sym_cls = "ok" if det > 0 else "fail"
        cmp_rows_html.append(
            f"<tr>"
            f"{_td(row['expected_name'])}"
            f"{_td(row['detection_tier'])}"
            f"{_td(_pct(row['expected_fraction']))}"
            f"{_td(_pct(row.get('observed_fraction', 0) or 0))}"
            f"{_td(_fc(row.get('log2_fold_change')))}"
            f"{_td(row.get('best_tier', '—'))}"
            f"{_td(sym, sym_cls)}"
            f"</tr>"
        )

    # Build unexpected-organism rows
    unexp_rows_html = []
    for row in top_unexpected.iter_rows(named=True):
        unexp_rows_html.append(
            f"<tr>"
            f"{_td(row['group_name'])}"
            f"{_td(row.get('kingdom', ''))}"
            f"{_td(row['n_reads_total'])}"
            f"{_td(row['best_tier'])}"
            f"{_td(row.get('evidence_calls', ''))}"
            f"</tr>"
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TPOS QC — Run {run_id}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         background: #f8fafc; color: #1e293b; margin: 0; padding: 24px; }}
  h1   {{ font-size: 1.4rem; font-weight: 700; margin-bottom: 4px; }}
  h2   {{ font-size: 1.05rem; font-weight: 600; margin: 24px 0 8px; color: #475569; }}
  .meta {{ color: #64748b; font-size: 0.85rem; margin-bottom: 24px; }}
  .badge {{ display: inline-block; padding: 4px 14px; border-radius: 999px;
            font-weight: 700; font-size: 1rem; color: #fff;
            background: {flag_color}; margin-bottom: 20px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.875rem; }}
  th, td {{ padding: 6px 10px; border: 1px solid #e2e8f0; text-align: left; }}
  th {{ background: #f1f5f9; font-weight: 600; }}
  tr:nth-child(even) {{ background: #f8fafc; }}
  .ok   {{ color: #16a34a; font-weight: 700; text-align: center; }}
  .fail {{ color: #dc2626; font-weight: 700; text-align: center; }}
  .summary-grid {{ display: grid; grid-template-columns: repeat(4, 1fr);
                   gap: 12px; margin-bottom: 24px; }}
  .card {{ background: #fff; border: 1px solid #e2e8f0; border-radius: 8px;
           padding: 12px 16px; }}
  .card-label {{ font-size: 0.75rem; color: #64748b; margin-bottom: 4px; }}
  .card-value {{ font-size: 1.15rem; font-weight: 700; }}
</style>
</head>
<body>
<h1>TPOS QC Report &mdash; ZymoBIOMICS Standard II Log</h1>
<div class="meta">Run {run_id} &nbsp;|&nbsp; {today} &nbsp;|&nbsp; lib_type: {r['lib_type']}</div>
<div class="badge">{qc_flag}</div>

<div class="summary-grid">
  <div class="card">
    <div class="card-label">MUST_DETECT rate</div>
    <div class="card-value">{must_rate:.0%}</div>
  </div>
  <div class="card">
    <div class="card-label">Spearman &rho; (log)</div>
    <div class="card-value">{f"{rho:.3f}" if rho is not None else "n/a"}</div>
  </div>
  <div class="card">
    <div class="card-label">LOD (lowest detected)</div>
    <div class="card-value">{f"{lod:.1e}" if lod is not None else "nothing"}</div>
  </div>
  <div class="card">
    <div class="card-label">Organisms detected</div>
    <div class="card-value">{n_det} / {n_exp}</div>
  </div>
</div>

<h2>Expected organism results (high &#x2192; low abundance)</h2>
<table>
  <thead>
    <tr>
      <th>Organism</th><th>Tier</th><th>Expected</th>
      <th>Observed</th><th>log2 FC</th><th>Best tier</th><th>Det.</th>
    </tr>
  </thead>
  <tbody>
    {"".join(cmp_rows_html)}
  </tbody>
</table>

<h2>Top unexpected organisms (not in ZymoBIOMICS panel, top 10 by read count)</h2>
{"<p><em>None detected.</em></p>" if not unexp_rows_html else f'''<table>
  <thead>
    <tr><th>Organism</th><th>Kingdom</th><th>Reads</th><th>Best tier</th><th>Evidence</th></tr>
  </thead>
  <tbody>
    {"".join(unexp_rows_html)}
  </tbody>
</table>'''}

<p style="margin-top:32px; color:#94a3b8; font-size:0.75rem;">
  Generated by analyse_report.py | CUB pipeline | CBAM / CHRU Brest
</p>
</body>
</html>"""

    Path(out_path).write_text(html, encoding="utf-8")
