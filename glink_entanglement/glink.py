import argparse
import os
import tempfile
import time
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "xdg-cache"))

import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from topoly import lasso_type


def point_rounding(num: float, threshold: float = 0.6) -> int:
    """Round absolute GLN values using the same threshold rule as the GE code."""
    if len(str(num).split("e")) == 2:
        num = 0.0

    if num % 1 >= threshold:
        return round(num)

    return int(str(num).split(".")[0])


def build_glink_matrix(ca_positions: np.ndarray) -> np.ndarray:
    """
    Calculate the Gaussian-linking dot matrix using GLink_ref.py broadcasting.

    Rows and columns correspond to C-alpha chain segments, so the matrix shape is
    (n_residues - 1, n_residues - 1).
    """
    ca_positions = ca_positions.astype(np.float32, copy=False)
    R = 0.5 * (ca_positions[:-1] + ca_positions[1:])
    dR = ca_positions[1:] - ca_positions[:-1]

    dR_cross = np.cross(dR[:, np.newaxis, :], dR[np.newaxis, :, :])
    diff = R[:, np.newaxis, :] - R[np.newaxis, :, :]
    norm_sq = np.einsum("ijk,ijk->ij", diff, diff)

    with np.errstate(divide="ignore", invalid="ignore"):
        inv_norm_cubed = np.where(norm_sq > 0.0, 1.0 / (norm_sq * np.sqrt(norm_sq)), 0.0)

    dot_matrix = np.einsum("ijk,ijk,ij->ij", diff, dR_cross, inv_norm_cubed)
    dot_matrix = np.nan_to_num(dot_matrix, copy=False)

    return dot_matrix / (4.0 * np.pi)


def build_prefix_matrix(matrix: np.ndarray) -> np.ndarray:
    """Build a padded 2D prefix-sum matrix for fast rectangular window sums."""
    return np.pad(matrix.cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)))


def window_sum(prefix_matrix: np.ndarray, row_start: int, row_stop: int, col_start: int, col_stop: int) -> float:
    """Return sum(matrix[row_start:row_stop, col_start:col_stop])."""
    return float(
        prefix_matrix[row_stop, col_stop]
        - prefix_matrix[row_start, col_stop]
        - prefix_matrix[row_stop, col_start]
        + prefix_matrix[row_start, col_start]
    )


def residue_key(residue) -> tuple:
    """Stable key for one residue inside a PDB chain/segment."""
    segid = getattr(residue, "segid", "") or ""
    chain_id = getattr(residue, "chainID", "") or ""
    chain = segid.strip() or chain_id.strip() or "SYSTEM"
    return chain, int(residue.resid), residue.resname


def atom_chain_key(atom) -> str:
    segid = getattr(atom, "segid", "") or ""
    chain_id = getattr(atom, "chainID", "") or ""
    return segid.strip() or chain_id.strip() or "SYSTEM"


def find_heavy_atom_contacts(chain_atoms, ca_atoms, cutoff: float = 4.5, min_sequence_separation: int = 4) -> list:
    """
    Return residue-index native contacts for one chain.

    A contact is present when any heavy-atom pair from two residues is within
    cutoff Angstrom and residue indices are at least min_sequence_separation
    positions apart.
    """
    residues = list(ca_atoms.residues)
    resid_to_index = {residue.resid: idx for idx, residue in enumerate(residues)}

    mapped_positions = []
    mapped_residue_indices = []
    for atom in chain_atoms:
        if atom.resid not in resid_to_index:
            continue
        mapped_positions.append(atom.position)
        mapped_residue_indices.append(resid_to_index[atom.resid])

    if not mapped_positions:
        return []

    mapped_positions = np.asarray(mapped_positions, dtype=np.float32)
    atom_residue_indices = np.asarray(mapped_residue_indices, dtype=np.int64)
    contact_pairs = set()

    for atom_i, atom_j in cKDTree(mapped_positions).query_pairs(cutoff):
        res_i = int(atom_residue_indices[atom_i])
        res_j = int(atom_residue_indices[atom_j])

        if res_i == res_j:
            continue

        if res_i > res_j:
            res_i, res_j = res_j, res_i

        if res_j - res_i >= min_sequence_separation:
            contact_pairs.add((res_i, res_j))

    return sorted(contact_pairs)


def map_topoly_crossings(crossings: list, residues: list) -> str:
    """
    Convert Topoly crossing labels from structure indices to PDB residue labels.

    Topoly documentation examples use 1-based structure indices, so `+25` maps
    to `residues[24]`. The chirality sign is preserved.
    """
    mapped_crossings = []
    for crossing in crossings:
        if not crossing:
            continue
        chirality = crossing[0]
        residue_index = int(crossing[1:]) - 1
        if residue_index < 0 or residue_index >= len(residues):
            continue
        _, resid, _ = residue_key(residues[residue_index])
        mapped_crossings.append(f"{chirality}{resid}")

    return " ".join(mapped_crossings)


def add_topoly_crossings(
    results: list,
    ca_positions: np.ndarray,
    residues: list,
    density: int = 0,
    min_dist: tuple = (10, 6, 5),
) -> list:
    """Add Topoly crossing residue columns and keep only contacts with crossings."""
    if not results:
        return []

    loop_indices = [
        [result["contact_i_index"] + 1, result["contact_j_index"] + 1]
        for result in results
    ]
    topoly_data = lasso_type(
        ca_positions.tolist(),
        loop_indices=loop_indices,
        more_info=True,
        density=density,
        min_dist=list(min_dist),
    )

    filtered_results = []
    for result in results:
        loop_key = (result["contact_i_index"] + 1, result["contact_j_index"] + 1)
        loop_data = topoly_data.get(loop_key, {})
        if not isinstance(loop_data, dict):
            loop_data = {}
        crossings_n = loop_data.get("crossingsN") or []
        crossings_c = loop_data.get("crossingsC") or []

        if result["Gn"] == 0:
            crossings_n = []
        if result["Gc"] == 0:
            crossings_c = []

        mapped_n = map_topoly_crossings(crossings_n, residues)
        mapped_c = map_topoly_crossings(crossings_c, residues)
        if not mapped_n and not mapped_c:
            continue

        result["crossingsN"] = mapped_n
        result["crossingsC"] = mapped_c
        result["crossings"] = " ".join(value for value in (mapped_n, mapped_c) if value)
        filtered_results.append(result)

    return filtered_results


def calculate_chain_glink(
    chain_atoms,
    cutoff: float = 4.5,
    threshold: float = 0.6,
    topoly_density: int = 0,
    topoly_min_dist: tuple = (10, 6, 5),
    nterm_threshold: int = 5,
    cterm_threshold: int = 5,
) -> list:
    """
    Calculate Gaussian linking values for heavy-atom contacts in one chain.

    The chain is represented by a heavy-atom AtomGroup. C-alpha atoms define the
    residue path used for the GLN calculation, while all heavy atoms are used to
    identify residue contacts. A contact is any residue pair whose heavy atoms
    come within `cutoff` Angstrom and whose sequence separation is at least 4.

    For each contact `(i, j)`, this function treats residues `i..j` as the loop
    and calculates two partial Gaussian linking values: `gn` against the
    N-terminal side of the chain and `gc` against the C-terminal side. Contacts
    with nonzero rounded `Gn` or `Gc` are Topoly candidates, but they are only
    retained in the final output if Topoly finds at least one crossing on a side
    whose rounded GLN value is nonzero.
    """
    ca_atoms = chain_atoms.select_atoms("name CA")
    if len(ca_atoms) < 2:
        return []

    contacts = find_heavy_atom_contacts(chain_atoms, ca_atoms, cutoff=cutoff)
    if not contacts:
        return []

    residues = list(ca_atoms.residues)
    glink_prefix = build_prefix_matrix(build_glink_matrix(ca_atoms.positions))
    n_residues = len(ca_atoms)
    results = []

    for i, j in contacts:
        if i - 5 > nterm_threshold:
            gn = window_sum(glink_prefix, nterm_threshold, i - 5, i, j)
        else:
            gn = 0.0

        cterm_start = j + 6
        cterm_stop = n_residues - (cterm_threshold + 1)
        if cterm_stop > cterm_start:
            gc = window_sum(glink_prefix, i, j, cterm_start, cterm_stop)
        else:
            gc = 0.0

        chain_i, resid_i, resname_i = residue_key(residues[i])
        chain_j, resid_j, resname_j = residue_key(residues[j])

        rounded_gn = point_rounding(abs(gn), threshold)
        rounded_gc = point_rounding(abs(gc), threshold)

        if rounded_gn == 0 and rounded_gc == 0:
            continue

        results.append(
            {
                "chain": chain_i,
                "resid_i": resid_i,
                "resname_i": resname_i,
                "resid_j": resid_j,
                "resname_j": resname_j,
                "contact_i_index": i,
                "contact_j_index": j,
                "gn": gn,
                "gc": gc,
                "Gn": rounded_gn,
                "Gc": rounded_gc,
            }
        )

    return add_topoly_crossings(
        results,
        ca_atoms.positions,
        residues,
        density=topoly_density,
        min_dist=topoly_min_dist,
    )


def calculate_pdb_glink(
    pdb_file: str,
    threshold: float = 0.6,
    cutoff: float = 4.5,
    topoly_density: int = 0,
    topoly_min_dist: tuple = (10, 6, 5),
) -> pd.DataFrame:
    """
    Calculate Gaussian linking values for heavy-atom contacts in a PDB file.

    This is the PDB-level driver. It loads the structure, selects protein heavy
    atoms with `protein and not name H*`, groups atoms by chain/segment, and
    calls `calculate_chain_glink()` for each chain.

    GLN values are calculated chain-by-chain because loop, N-terminal, and
    C-terminal windows are sequence-based. The returned DataFrame contains only
    contacts with both nonzero rounded `Gn` or `Gc` and at least one matching
    Topoly crossing residue, with residue IDs, residue names, zero-based contact
    indices, raw `gn`/`gc`, rounded `Gn`/`Gc`, and Topoly crossing residues.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"1 A\^3 CRYST1 record.*",
            category=UserWarning,
            module=r"MDAnalysis\.coordinates\.PDB",
        )
        universe = mda.Universe(pdb_file, format="PDB")
    heavy_atoms = universe.select_atoms("protein and not name H*")
    all_results = []

    chain_to_atom_indices = {}
    for atom in heavy_atoms:
        chain_to_atom_indices.setdefault(atom_chain_key(atom), []).append(atom.index)

    for chain in sorted(chain_to_atom_indices):
        chain_atoms = universe.atoms[chain_to_atom_indices[chain]]
        all_results.extend(
            calculate_chain_glink(
                chain_atoms,
                cutoff=cutoff,
                threshold=threshold,
                topoly_density=topoly_density,
                topoly_min_dist=topoly_min_dist,
            )
        )

    columns = [
        "chain",
        "resid_i",
        "resname_i",
        "resid_j",
        "resname_j",
        "contact_i_index",
        "contact_j_index",
        "gn",
        "gc",
        "Gn",
        "Gc",
        "crossingsN",
        "crossingsC",
        "crossings",
    ]
    return pd.DataFrame(all_results, columns=columns)


def resolve_output_path(pdb_file: str, output: str | None) -> Path:
    default_name = f"{Path(pdb_file).stem}_glink_contacts.csv"

    if output is None:
        return Path(default_name)

    output_path = Path(output)
    if output_path.suffix.lower() == ".csv":
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return output_path

    output_path.mkdir(parents=True, exist_ok=True)
    return output_path / default_name


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate per-contact Gaussian linking values from an all-atom PDB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--PDB", required=True, help="Required. Input all-atom PDB file.")
    parser.add_argument(
        "-o",
        "--output",
        help="Optional. Output directory or CSV path. Directory writes <pdb>_glink_contacts.csv inside it.",
    )
    parser.add_argument(
        "--GLN_threshold",
        type=float,
        default=0.6,
        help="Optional. Threshold used to round absolute gn/gc values into Gn/Gc.",
    )
    parser.add_argument(
        "--contact_cutoff",
        type=float,
        default=4.5,
        help="Optional. Heavy-atom distance cutoff in Angstrom for residue contacts.",
    )
    parser.add_argument(
        "--topoly_density",
        type=int,
        default=0,
        help="Optional. Topoly minimal-surface triangulation density. Use 1 for Topoly's default precision.",
    )
    parser.add_argument(
        "--topoly_min_dist",
        type=int,
        nargs=3,
        default=(10, 6, 5),
        metavar=("CROSSING", "LOOP", "TAIL_END"),
        help="Optional. Topoly crossing-reduction distances.",
    )
    args = parser.parse_args()

    output = resolve_output_path(args.PDB, args.output)

    total_start = time.perf_counter()
    data = calculate_pdb_glink(
        args.PDB,
        threshold=args.GLN_threshold,
        cutoff=args.contact_cutoff,
        topoly_density=args.topoly_density,
        topoly_min_dist=tuple(args.topoly_min_dist),
    )
    output_data = data.copy()
    output_data.insert(
        0,
        "contact",
        output_data["resname_i"]
        + output_data["resid_i"].astype(str)
        + "-"
        + output_data["resname_j"]
        + output_data["resid_j"].astype(str),
    )
    output_data = output_data.drop(
        columns=["chain", "resid_i", "resname_i", "resid_j", "resname_j"]
    )
    output_data = output_data.rename(
        columns={"contact_i_index": "i", "contact_j_index": "j"}
    )
    for column in ("gn", "gc"):
        output_data[column] = output_data[column].map(lambda value: f"{value:.3f}")
    output_data.to_csv(output, index=False)
    elapsed = time.perf_counter() - total_start

    print(f"Wrote {len(data)} contacts to {output}")
    print(f"Runtime: {elapsed:.6f} seconds")


if __name__ == "__main__":
    main()
