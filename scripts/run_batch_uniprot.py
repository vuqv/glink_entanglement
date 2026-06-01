#!/usr/bin/env python3
"""Batch wrapper for a table with a Uniprot column.

Each Uniprot value is resolved to <pdb_dir>/<Uniprot>.pdb, then processed with
one PDB per worker. Raw GLINK output is written first; clustering runs only if
raw entanglements are present.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import time

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

def infer_delimiter(table_file: Path, explicit_delimiter: str | None) -> str:
    if explicit_delimiter is not None:
        return explicit_delimiter
    if table_file.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    return ","


def unique_nonempty_values(values) -> list[str]:
    uniprots = []
    seen = set()
    for value in values:
        if pd.isna(value):
            continue
        value = str(value).strip()
        if not value or value in seen:
            continue
        seen.add(value)
        uniprots.append(value)
    return uniprots


def read_uniprots(table_file: Path, column: str = "Uniprot", delimiter: str | None = None) -> list[str]:
    if table_file.suffix.lower() in {".pkl", ".pickle"}:
        data = pd.read_pickle(table_file)
        if column not in data.columns:
            raise ValueError(f"Column {column!r} not found in {table_file}. Found columns: {list(data.columns)}")
        return unique_nonempty_values(data[column])

    delimiter = infer_delimiter(table_file, delimiter)
    with table_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {table_file}")
        if column not in reader.fieldnames:
            raise ValueError(f"Column {column!r} not found in {table_file}. Found columns: {reader.fieldnames}")

        return unique_nonempty_values(row.get(column) for row in reader)


def resolve_uniprot_pdbs(uniprots: list[str], pdb_dir: Path, extension: str = ".pdb") -> tuple[list[Path], list[dict]]:
    if not extension.startswith("."):
        extension = f".{extension}"

    pdb_files = []
    missing_rows = []
    for uniprot in uniprots:
        pdb_file = pdb_dir / f"{uniprot}{extension}"
        if pdb_file.exists():
            pdb_files.append(pdb_file)
        else:
            missing_rows.append(
                {
                    "pdb": str(pdb_file),
                    "raw_csv": "",
                    "clustered_csv": "",
                    "status": "missing_pdb",
                    "clustered": False,
                    "seconds": "0.000",
                    "error": f"PDB file not found for Uniprot {uniprot}",
                }
            )

    return pdb_files, missing_rows


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
        description="Run glink batch processing from a table containing a Uniprot column.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/run_batch_uniprot.py -t proteins.pkl -p PDB -o batch_results -j 32\n"
            "  python scripts/run_batch_uniprot.py -t proteins.csv -p PDB -o batch_results -j 32\n"
            "  python scripts/run_batch_uniprot.py -t proteins.tsv -p PDB -o batch_results --delimiter '\\t'\n"
            "  python scripts/run_batch_uniprot.py -t proteins.pkl -p PDB -o batch_results --column UniProtID"
        ),
    )
    parser.add_argument("-t", "--table", required=True, help="Required. CSV/TSV or pickle DataFrame containing a Uniprot column.")
    parser.add_argument("-p", "--pdb_dir", required=True, help="Required. Directory containing <Uniprot>.pdb files.")
    parser.add_argument("-o", "--outdir", required=True, help="Required. Batch output directory.")
    parser.add_argument("--column", default="Uniprot", help="Column containing Uniprot IDs.")
    parser.add_argument("--delimiter", help="Input table delimiter. Defaults to comma, or tab for .tsv/.tab files.")
    parser.add_argument("--extension", default=".pdb", help="PDB filename extension appended to each Uniprot ID.")
    parser.add_argument("-j", "--workers", type=int, default=os.cpu_count() or 1, help="Number of worker processes.")
    parser.add_argument("--force", action="store_true", help="Recompute outputs even if CSV files already exist.")
    parser.add_argument("--summary", help="Summary CSV path. Defaults to <outdir>/batch_summary.csv.")
    args = parser.parse_args()

    table_file = Path(args.table)
    pdb_dir = Path(args.pdb_dir)
    outdir = Path(args.outdir)
    raw_dir = outdir / "raw"
    clustered_dir = outdir / "clustered"
    raw_dir.mkdir(parents=True, exist_ok=True)
    clustered_dir.mkdir(parents=True, exist_ok=True)

    uniprots = read_uniprots(table_file, column=args.column, delimiter=args.delimiter)
    if not uniprots:
        raise SystemExit(f"No Uniprot values found in {table_file}")

    pdb_files, rows = resolve_uniprot_pdbs(uniprots, pdb_dir, extension=args.extension)
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
    if work_items:
        max_workers = min(args.workers, len(work_items))
        print(f"Processing {len(work_items)} PDB files with {max_workers} workers")
        if rows:
            print(f"Skipping {len(rows)} missing PDB files")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_one_pdb, item): item[0] for item in work_items}
            for future in as_completed(futures):
                result = future.result()
                rows.append(result)
                print(f"[{result['status']}] {result['pdb']} ({result['seconds']} s)")
    else:
        print("No existing PDB files to process")

    rows.sort(key=lambda row: row["pdb"])
    write_summary(summary_file, rows)
    elapsed = time.perf_counter() - start
    failed = sum(row["status"] == "failed" for row in rows)
    missing = sum(row["status"] == "missing_pdb" for row in rows)
    clustered = sum(row["clustered"] for row in rows)
    print(f"Summary: {summary_file}")
    print(f"Clustered: {clustered}; missing: {missing}; failed: {failed}; runtime: {elapsed:.3f} seconds")

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
