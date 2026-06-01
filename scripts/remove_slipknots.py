#!/usr/bin/env python3
"""Remove slipknot crossing pairs from clustered GLINK output.

The input is a clustered CSV from glink-cluster. N-terminal crossings are
checked from larger residue number to smaller residue number; C-terminal
crossings are checked from smaller residue number to larger residue number.
Adjacent crossings with opposite signs cancel, and cancellation is repeated until
no adjacent opposite-sign pair remains.
"""

from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

import pandas as pd


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


CROSSING_PATTERN = re.compile(r"^[+-]\d+$")


def parse_crossings(value) -> list[str]:
    if pd.isna(value):
        return []

    crossings = []
    for token in str(value).replace(",", " ").split():
        token = token.strip()
        if CROSSING_PATTERN.fullmatch(token):
            crossings.append(token)
    return crossings


def crossing_number(crossing: str) -> int:
    return int(crossing[1:])


def crossing_sign(crossing: str) -> str:
    return crossing[0]


def order_crossings(crossings: list[str], side: str) -> list[str]:
    reverse = side.upper() == "N"
    return sorted(crossings, key=crossing_number, reverse=reverse)


def remove_canceling_pairs(ordered_crossings: list[str]) -> list[str]:
    stack: list[str] = []
    for crossing in ordered_crossings:
        if stack and crossing_sign(stack[-1]) != crossing_sign(crossing):
            stack.pop()
        else:
            stack.append(crossing)
    return stack


def remove_slipknots_from_crossings(value, side: str) -> str:
    crossings = parse_crossings(value)
    ordered = order_crossings(crossings, side)
    reduced = remove_canceling_pairs(ordered)
    return " ".join(reduced)


def resolve_output_path(input_csv: Path, output: str | None) -> Path:
    default_name = f"{input_csv.stem}_no_slipknot.csv"
    if output is None:
        return input_csv.with_name(default_name)

    output_path = Path(output)
    if output_path.suffix.lower() == ".csv":
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return output_path

    output_path.mkdir(parents=True, exist_ok=True)
    return output_path / default_name


def remove_slipknots(input_csv: str, output: str | None = None, keep_empty: bool = False) -> pd.DataFrame:
    input_path = Path(input_csv)
    data = pd.read_csv(
        input_path,
        keep_default_na=False,
        dtype={"crossingsN": str, "crossingsC": str, "crossings": str},
    )
    required_columns = {"crossingsN", "crossingsC", "crossings"}
    missing = required_columns.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns in {input_csv}: {sorted(missing)}")

    processed = data.copy()
    processed = processed[~processed.astype(str).apply(lambda row: all(value.strip() == "" for value in row), axis=1)].copy()
    processed["crossingsN"] = processed["crossingsN"].map(lambda value: remove_slipknots_from_crossings(value, "N"))
    processed["crossingsC"] = processed["crossingsC"].map(lambda value: remove_slipknots_from_crossings(value, "C"))
    processed["crossings"] = processed.apply(
        lambda row: " ".join(value for value in (row["crossingsN"], row["crossingsC"]) if value),
        axis=1,
    )
    if not keep_empty:
        processed = processed[processed["crossings"] != ""].copy()

    output_path = resolve_output_path(input_path, output)
    processed.to_csv(output_path, index=False)
    return processed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove slipknot crossing pairs from clustered GLINK CSV output.",
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  python scripts/remove_slipknots.py -i clustered/2ww4_glink_contacts_clustered.csv\n"
            "  python scripts/remove_slipknots.py -i clustered.csv -o cleaned\n"
            "  python scripts/remove_slipknots.py -i clustered.csv -o cleaned/no_slipknot.csv"
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="Required. Clustered GLINK CSV with crossingsN/crossingsC columns.")
    parser.add_argument(
        "-o",
        "--output",
        help="Optional. Output directory or CSV path. Defaults to <input_stem>_no_slipknot.csv beside input.",
    )
    parser.add_argument(
        "--keep_empty",
        action="store_true",
        help="Keep rows whose crossings all cancel. By default these rows are dropped.",
    )
    args = parser.parse_args()

    start = time.perf_counter()
    output_path = resolve_output_path(Path(args.input), args.output)
    processed = remove_slipknots(args.input, str(output_path), keep_empty=args.keep_empty)
    elapsed = time.perf_counter() - start

    print(f"Wrote {len(processed)} rows to {output_path}")
    print(f"Runtime: {elapsed:.6f} seconds")


if __name__ == "__main__":
    main()
