#!/usr/bin/env python3
"""Batch wrapper for glink-entanglement.

Each worker process handles one PDB file:
1. run raw GLINK entanglement calculation;
2. if the raw CSV has at least one result row, run clustering.

This script is intentionally separate from the package CLI because batch runs are
cluster/job-scheduler specific and are not always needed.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


PDB_SUFFIXES = {".pdb", ".ent", ".cif"}


def discover_pdb_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]

    pdb_files = [
        path
        for path in input_path.iterdir()
        if path.is_file() and path.suffix.lower() in PDB_SUFFIXES
    ]
    return sorted(pdb_files)


def csv_has_rows(csv_file: Path) -> bool:
    if not csv_file.exists() or csv_file.stat().st_size == 0:
        return False

    with csv_file.open(newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        return any(True for _ in reader)


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


def process_one_pdb(args_tuple: tuple) -> dict:
    (
        pdb_file,
        raw_dir,
        clustered_dir,
        force,
    ) = args_tuple

    pdb_file = Path(pdb_file)
    raw_dir = Path(raw_dir)
    clustered_dir = Path(clustered_dir)
    raw_csv = raw_dir / f"{pdb_file.stem}_glink_contacts.csv"
    clustered_csv = clustered_dir / f"{raw_csv.stem}_clustered.csv"

    start = time.perf_counter()
    status = "ok"
    clustered = False
    error = ""

    try:
        if force or not raw_csv.exists():
            glink_command = [
                sys.executable,
                "-m",
                "glink_entanglement",
                "-f",
                str(pdb_file),
                "-o",
                str(raw_csv),
            ]
            run_command(glink_command)

        if csv_has_rows(raw_csv):
            if force or not clustered_csv.exists():
                cluster_command = [
                    sys.executable,
                    "-m",
                    "glink_entanglement.clustering",
                    "-f",
                    str(raw_csv),
                    "-o",
                    str(clustered_dir),
                    "-w",
                    str(clustered_csv),
                ]
                run_command(cluster_command)
            clustered = True
        else:
            status = "no_entanglement"
    except subprocess.CalledProcessError as exc:
        status = "failed"
        error = f"command failed with exit code {exc.returncode}: {' '.join(exc.cmd)}"
    except Exception as exc:  # noqa: BLE001 - batch summaries should capture per-PDB failures.
        status = "failed"
        error = repr(exc)

    return {
        "pdb": str(pdb_file),
        "raw_csv": str(raw_csv),
        "clustered_csv": str(clustered_csv) if clustered else "",
        "status": status,
        "clustered": clustered,
        "seconds": f"{time.perf_counter() - start:.3f}",
        "error": error,
    }


def write_summary(summary_file: Path, rows: list[dict]) -> None:
    fieldnames = ["pdb", "raw_csv", "clustered_csv", "status", "clustered", "seconds", "error"]
    with summary_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run glink and clustering over many PDB files, one PDB per worker.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/run_batch.py -i PDB -o batch_results -j 32\n"
            "  python scripts/run_batch.py -i one_structure.pdb -o batch_results --force"
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="Required. PDB file or directory of PDB files.")
    parser.add_argument("-o", "--outdir", required=True, help="Required. Batch output directory.")
    parser.add_argument("-j", "--workers", type=int, default=os.cpu_count() or 1, help="Number of worker processes.")
    parser.add_argument("--force", action="store_true", help="Recompute outputs even if CSV files already exist.")
    parser.add_argument("--summary", help="Summary CSV path. Defaults to <outdir>/batch_summary.csv.")
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    raw_dir = outdir / "raw"
    clustered_dir = outdir / "clustered"
    raw_dir.mkdir(parents=True, exist_ok=True)
    clustered_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = discover_pdb_files(input_path)
    if not pdb_files:
        raise SystemExit(f"No PDB files found in {input_path}")

    summary_file = Path(args.summary) if args.summary else outdir / "batch_summary.csv"
    work_items = [
        (
            pdb_file,
            raw_dir,
            clustered_dir,
            args.force,
        )
        for pdb_file in pdb_files
    ]

    start = time.perf_counter()
    rows = []
    max_workers = min(args.workers, len(work_items))
    print(f"Processing {len(work_items)} PDB files with {max_workers} workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_one_pdb, item): item[0] for item in work_items}
        for future in as_completed(futures):
            result = future.result()
            rows.append(result)
            print(f"[{result['status']}] {result['pdb']} ({result['seconds']} s)")

    rows.sort(key=lambda row: row["pdb"])
    write_summary(summary_file, rows)
    elapsed = time.perf_counter() - start
    failed = sum(row["status"] == "failed" for row in rows)
    clustered = sum(row["clustered"] for row in rows)
    print(f"Summary: {summary_file}")
    print(f"Clustered: {clustered}; failed: {failed}; runtime: {elapsed:.3f} seconds")

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
