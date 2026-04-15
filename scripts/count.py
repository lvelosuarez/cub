#!/usr/bin/env python3
"""
Count reads in FASTQ(.gz) files recursively, using multithreading.

- By default: only counts "R1" / "1" files (first mate).
- Match examples like:
    CasClinique7_ADN_R1.fastq.gz
    CasClinique7_ADN_dedup_R1.fastq.gz
    CasClinique7_ADN_dedup_R1.clean_1.fastq.gz
    sample_1.fastq.gz, sample_1.fq, etc.

Output: TSV with columns:
    file    output
where "file" is the relative path from the input directory.
"""

import argparse
import gzip
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def is_fastq(path: Path) -> bool:
    """
    Return True if the file looks like a FASTQ(.gz).
    """
    name = path.name.lower()
    return name.endswith(".fastq") or name.endswith(".fq") or name.endswith(".fastq.gz") or name.endswith(".fq.gz")


def is_r1_like(path: Path) -> bool:
    """
    Return True if the file name looks like R1 / first mate.

    Examples matched:
    - *_R1.fastq.gz
    - *_R1_clean_1.fastq.gz
    - *_1.fastq.gz / *_1.fq.gz / *_1.fastq / *_1.fq
    """
    name = path.name
    if "_R1" in name:
        return True
    # generic Illumina-style *_1.fastq etc.
    if name.endswith(("_1.fastq", "_1.fq", "_1.fastq.gz", "_1.fq.gz")):
        return True
    return False


def count_fastq_reads(path: Path) -> int:
    """
    Count reads in a FASTQ/FASTQ.GZ file.
    FASTQ: 4 lines per read.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    # handle .fastq.gz with double suffix
    if path.name.endswith(".fastq.gz") or path.name.endswith(".fq.gz"):
        opener = gzip.open

    line_count = 0
    try:
        with opener(path, "rt") as f:
            for _ in f:
                line_count += 1
    except Exception as e:
        print(f"[ERROR] Failed to read {path}: {e}", file=sys.stderr)
        return -1

    return line_count // 4


def main():
    parser = argparse.ArgumentParser(description="Recursively count reads in FASTQ(.gz) files using multithreading.")
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Root directory to scan recursively.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file (file, output).",
    )
    parser.add_argument(
        "--all-fastq",
        action="store_true",
        help="If set, count ALL fastq files (R1 and R2). " "By default, only R1/1 files are counted.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 4,
        help="Number of worker threads (default: number of CPUs).",
    )

    args = parser.parse_args()
    root = Path(args.input_dir).resolve()

    if not root.is_dir():
        raise SystemExit(f"[ERROR] Not a directory: {root}")

    print(f"[INFO] Scanning {root} recursively for FASTQ files...")
    candidates = []

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if not is_fastq(p):
            continue
        if not args.all_fastq and not is_r1_like(p):
            continue
        candidates.append(p)

    if not candidates:
        print("[WARN] No matching FASTQ files found.", file=sys.stderr)
        # still write an empty file with header
        with open(args.output, "w") as out:
            out.write("file\toutput\n")
        return

    print(f"[INFO] Found {len(candidates)} files to count.")
    print(f"[INFO] Using {args.threads} threads.")

    results = []

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_path = {executor.submit(count_fastq_reads, p): p for p in candidates}

        for future in as_completed(future_to_path):
            p = future_to_path[future]
            rel = p.relative_to(root)
            try:
                reads = future.result()
            except Exception as e:
                print(f"[ERROR] Exception counting {p}: {e}", file=sys.stderr)
                reads = -1
            results.append((str(rel), reads))

    # sort results for nice output: by path
    results.sort(key=lambda x: x[0])

    print(f"[INFO] Writing output to {args.output}")
    with open(args.output, "w") as out:
        out.write("file\toutput\n")
        for path_str, reads in results:
            out.write(f"{path_str}\t{reads}\n")

    print("[DONE]")


if __name__ == "__main__":
    main()
