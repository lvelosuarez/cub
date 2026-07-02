#!/usr/bin/env python3
"""
generate_report.py

Builds the Patient Results Explorer HTML report for a single run.

Reads:
  - a single .parquet file containing all patient/sample rows for the run
  - (optional) a JSON file with a TPOS QC summary block

Writes:
  - a single self-contained HTML file (template + injected data)

Usage:
  python generate_report.py \
      --parquet /path/to/run14.parquet \
      --template /path/to/template.html \
      --output /path/to/output/run14_report.html \
      [--qc-summary /path/to/qc_summary.json]
"""

import argparse
import json
import sys
from pathlib import Path

import polars as pl

# Kingdom values the table filter cares about. Anything else is shown
# regardless of the kingdom filter state (see template behavior).
KNOWN_KINGDOMS = ["Bacteria", "Eukaryota", "Viruses"]


def parse_args():
    p = argparse.ArgumentParser(description="Generate the Patient Results Explorer HTML report.")
    p.add_argument("--parquet", required=True, type=Path,
                    help="Path to the run's .parquet file (sample_id, group_name, kingdom, "
                         "decontam_decision, Library_size, microbial_reads, ...).")
    p.add_argument("--template", required=True, type=Path,
                    help="Path to template.html.")
    p.add_argument("--output", required=True, type=Path,
                    help="Full path (including filename) of the HTML file to create.")
    p.add_argument("--qc-summary", required=False, type=Path, default=None,
                    help="Optional path to a JSON file with the TPOS QC summary block.")
    return p.parse_args()


def load_patient_data(parquet_path: Path):
    """
    Reads the parquet file and groups rows by sample_id, preserving
    first-appearance order. Returns a list (not a dict) so that
    numeric-looking sample_id values can't get silently reordered by
    JS object key semantics.
    """
    df = pl.read_parquet(parquet_path)

    if "sample_id" not in df.columns:
        sys.exit(f"ERROR: column 'sample_id' not found in {parquet_path}")

    df = df.with_columns(pl.col("sample_id").cast(pl.Utf8))

    # Preserve first-appearance order of sample_id.
    seen_order = []
    seen_set = set()
    for sid in df["sample_id"].to_list():
        if sid not in seen_set:
            seen_set.add(sid)
            seen_order.append(sid)

    rows_by_id = {sid: [] for sid in seen_order}
    for row in df.iter_rows(named=True):
        rows_by_id[row["sample_id"]].append(row)

    patient_data = [
        {"sample_id": sid, "rows": rows_by_id[sid]}
        for sid in seen_order
    ]
    return patient_data, df, seen_order


def load_qc_summary(qc_summary_path: Path):
    if qc_summary_path is None:
        return None
    if not qc_summary_path.exists():
        sys.exit(f"ERROR: --qc-summary file not found: {qc_summary_path}")
    try:
        with open(qc_summary_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        sys.exit(f"ERROR: --qc-summary is not valid JSON ({qc_summary_path}): {e}")
    return data


def build_qc_reads(df: pl.DataFrame, seen_order):
    """
    Aggregates, per sample_id (same order as PATIENT_DATA):
      - library_size       (assumed constant across the sample's rows; first row used)
      - microbial_reads    (same assumption)
      - kingdom_pass        counts of PASS rows per known kingdom
    """
    has_lib = "Library_size" in df.columns
    has_micro = "microbial_reads" in df.columns
    has_kingdom = "kingdom" in df.columns
    has_decision = "decontam_decision" in df.columns

    qc_reads = []
    for sid in seen_order:
        sub = df.filter(pl.col("sample_id") == sid)

        library_size = sub["Library_size"][0] if has_lib and sub.height > 0 else None
        microbial_reads = sub["microbial_reads"][0] if has_micro and sub.height > 0 else None

        kingdom_pass = {}
        kingdom_total = {}
        for k in KNOWN_KINGDOMS:
            if has_kingdom:
                k_rows = sub.filter(pl.col("kingdom") == k)
                kingdom_total[k] = k_rows.height
                kingdom_pass[k] = k_rows.filter(pl.col("decontam_decision") == "PASS").height if has_decision else 0
            else:
                kingdom_total[k] = 0
                kingdom_pass[k] = 0

        qc_reads.append({
            "sample_id": sid,
            "library_size": library_size,
            "microbial_reads": microbial_reads,
            "kingdom_pass": kingdom_pass,
            "kingdom_total": kingdom_total,
        })
       
    return qc_reads


def to_json(obj):
    """JSON-dump with safety escaping so the data can't break out of the <script> tag."""
    s = json.dumps(obj, ensure_ascii=False, default=str)
    return s.replace("</", "<\\/")


def generate_html(template_path: Path, output_path: Path,
                   patient_data, qc_summary, qc_reads):
    template = template_path.read_text(encoding="utf-8")

    html = template
    html = html.replace("__PATIENT_DATA_JSON__", to_json(patient_data))
    html = html.replace("__QC_SUMMARY_JSON__", to_json(qc_summary))
    html = html.replace("__QC_READS_JSON__", to_json(qc_reads))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"\u2713 HTML report written to {output_path}")


def main():
    args = parse_args()

    print("Loading parquet...")
    patient_data, df, seen_order = load_patient_data(args.parquet)
    print(f"  {len(patient_data)} samples, {df.height} total rows.")

    print("Loading QC summary..." if args.qc_summary else "No --qc-summary provided, skipping.")
    qc_summary = load_qc_summary(args.qc_summary)

    print("Building per-sample reads table...")
    qc_reads = build_qc_reads(df, seen_order)

    print("Generating HTML...")
    generate_html(args.template, args.output, patient_data, qc_summary, qc_reads)


if __name__ == "__main__":
    main()
