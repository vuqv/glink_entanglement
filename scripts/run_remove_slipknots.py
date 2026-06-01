#!/usr/bin/env python3
"""Run slipknot removal over a directory of clustered GLINK CSV files."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import time
from pathlib import Path


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def load_remove_slipknots_function():
    script_path = Path(__file__).resolve().with_name("remove_slipknots.py")
    spec = importlib.util.spec_from_file_location("remove_slipknots", script_path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"Could not load {script_path}")
    spec.loader.exec_module(module)
    return module.remove_slipknots


def discover_clustered_csvs(input_dir: Path) -> list[Path]:
    return sorted(path for path in input_dir.glob("*.csv") if path.is_file())


def output_path_for(input_csv: Path, output_dir: Path) -> Path:
    return output_dir / f"{input_csv.stem}_no_slipknot.csv"


def write_summary(summary_file: Path, rows: list[dict]) -> None:
    fieldnames = ["input_csv", "output_csv", "status", "input_rows", "output_rows", "seconds", "error"]
    with summary_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Serially remove slipknot crossings from all clustered GLINK CSV files in a directory.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/run_remove_slipknots.py -i batch_results/clustered -o batch_results/no_slipknot\n"
            "  python scripts/run_remove_slipknots.py -b batch_results\n"
            "  python scripts/run_remove_slipknots.py -b batch_results --keep_empty"
        ),
    )
    parser.add_argument(
        "-b",
        "--batch_dir",
        default="batch_results",
        help="Batch result directory. Used to infer input/output dirs when -i/-o are not provided.",
    )
    parser.add_argument("-i", "--input_dir", help="Directory containing clustered CSV files. Defaults to <batch_dir>/clustered.")
    parser.add_argument("-o", "--output_dir", help="Output directory. Defaults to <batch_dir>/no_slipknot.")
    parser.add_argument("--keep_empty", action="store_true", help="Keep rows whose crossings all cancel. By default these rows are dropped.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing cleaned CSV files.")
    parser.add_argument("--summary", help="Summary CSV path. Defaults to <output_dir>/remove_slipknots_summary.csv.")
    args = parser.parse_args()

    batch_dir = Path(args.batch_dir)
    input_dir = Path(args.input_dir) if args.input_dir else batch_dir / "clustered"
    output_dir = Path(args.output_dir) if args.output_dir else batch_dir / "no_slipknot"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_file = Path(args.summary) if args.summary else output_dir / "remove_slipknots_summary.csv"

    csv_files = discover_clustered_csvs(input_dir)
    if not csv_files:
        raise SystemExit(f"No CSV files found in {input_dir}")

    remove_slipknots = load_remove_slipknots_function()
    rows = []
    total_start = time.perf_counter()
    print(f"Processing {len(csv_files)} clustered CSV files")

    for input_csv in csv_files:
        start = time.perf_counter()
        output_csv = output_path_for(input_csv, output_dir)
        status = "ok"
        error = ""
        input_rows = ""
        output_rows = ""

        try:
            if output_csv.exists() and not args.force:
                status = "skipped"
            else:
                processed = remove_slipknots(str(input_csv), str(output_csv), keep_empty=args.keep_empty)
                output_rows = len(processed)
                input_rows = sum(1 for _ in input_csv.open()) - 1
        except ValueError as exc:
            if "Missing required columns" in str(exc):
                status = "skipped_missing_columns"
                error = str(exc)
            else:
                status = "failed"
                error = repr(exc)
        except Exception as exc:  # noqa: BLE001 - summary should record per-file failures.
            status = "failed"
            error = repr(exc)

        seconds = time.perf_counter() - start
        rows.append(
            {
                "input_csv": str(input_csv),
                "output_csv": str(output_csv),
                "status": status,
                "input_rows": input_rows,
                "output_rows": output_rows,
                "seconds": f"{seconds:.3f}",
                "error": error,
            }
        )
        print(f"[{status}] {input_csv} -> {output_csv} ({seconds:.3f} s)")

    write_summary(summary_file, rows)
    elapsed = time.perf_counter() - total_start
    failed = sum(row["status"] == "failed" for row in rows)
    skipped = sum(row["status"] in {"skipped", "skipped_missing_columns"} for row in rows)
    print(f"Summary: {summary_file}")
    print(f"Failed: {failed}; skipped: {skipped}; runtime: {elapsed:.3f} seconds")

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
