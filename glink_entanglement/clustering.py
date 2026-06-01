import argparse
import time
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


ORGANISM_CUTOFFS = {
    "Human": 52,
    "Ecoli": 57,
    "Yeast": 49,
}


@dataclass(frozen=True)
class Entanglement:
    chain: str
    i: int
    j: int
    crossings: tuple
    gn: float
    gc: float
    Gn: int
    Gc: int
    num_contacts: int
    contacts: tuple


def crossing_number(crossing: str) -> int:
    return int(crossing[1:])


def sort_crossings(crossings: list[str]) -> tuple:
    return tuple(sorted(crossings, key=crossing_number))


def parse_crossings(value) -> tuple:
    if pd.isna(value):
        return ()
    crossings = str(value).replace(",", " ").split()
    crossings = [crossing.strip() for crossing in crossings if crossing.strip()]
    return sort_crossings(crossings)


def contact_label(i: int, j: int) -> str:
    return f"{i}-{j}"


def loops_overlap(ent1: Entanglement, ent2: Entanglement) -> bool:
    return (
        ent2.i <= ent1.i <= ent2.j
        or ent2.i <= ent1.j <= ent2.j
        or ent1.i <= ent2.i <= ent1.j
        or ent1.i <= ent2.j <= ent1.j
    )


def entanglement_coord(ent: Entanglement) -> tuple:
    return (ent.i, ent.j, *[crossing_number(crossing) for crossing in ent.crossings])


def same_chirality(ent1: Entanglement, ent2: Entanglement) -> bool:
    return [cr[0] for cr in ent1.crossings] == [cr[0] for cr in ent2.crossings]


def merge_into_representative(rep: Entanglement, merged: Entanglement) -> Entanglement:
    contacts = tuple(dict.fromkeys((*rep.contacts, *merged.contacts)))
    return replace(rep, num_contacts=rep.num_contacts + merged.num_contacts, contacts=contacts)


def read_glink_csv(csv_file: str) -> list[Entanglement]:
    """
    Read glink.py output and return rows with Topoly crossings.

    The standalone GLN calculator already filters contacts to nonzero rounded
    Gn/Gc. This clustering stage additionally requires at least one crossing
    residue, because representative entanglements are defined by crossing sets.
    """
    data = pd.read_csv(csv_file, dtype={"chain": str, "crossings": str}, keep_default_na=False)
    required_columns = {"chain", "resid_i", "resid_j", "gn", "gc", "Gn", "Gc", "crossings"}
    missing = required_columns.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns in {csv_file}: {sorted(missing)}")

    entanglements = []
    for row in data.itertuples(index=False):
        crossings = parse_crossings(getattr(row, "crossings"))
        if not crossings:
            continue

        i = int(getattr(row, "resid_i"))
        j = int(getattr(row, "resid_j"))
        if i > j:
            i, j = j, i

        entanglements.append(
            Entanglement(
                chain=str(getattr(row, "chain")),
                i=i,
                j=j,
                crossings=crossings,
                gn=float(getattr(row, "gn")),
                gc=float(getattr(row, "gc")),
                Gn=int(getattr(row, "Gn")),
                Gc=int(getattr(row, "Gc")),
                num_contacts=1,
                contacts=(contact_label(i, j),),
            )
        )

    return entanglements


def select_minimal_loop_by_crossing_set(entanglements: list[Entanglement]) -> list[Entanglement]:
    grouped = defaultdict(list)
    for ent in entanglements:
        grouped[(ent.chain, ent.crossings)].append(ent)

    representatives = []
    for group in grouped.values():
        rep = min(group, key=lambda ent: (ent.j - ent.i, ent.i, ent.j))
        contacts = tuple(contact_label(ent.i, ent.j) for ent in group)
        representatives.append(
            replace(rep, num_contacts=len(group), contacts=contacts)
        )

    return representatives


def merge_larger_overlapping_loops(entanglements: list[Entanglement]) -> list[Entanglement]:
    """
    Merge same-size, same-chirality crossing sets by keeping the smaller loop.

    This mirrors Step 3 in clustering.py.
    """
    active = list(entanglements)
    changed = True

    while changed:
        changed = False
        for idx_a in range(len(active)):
            if changed:
                break
            for idx_b in range(idx_a + 1, len(active)):
                ent_a = active[idx_a]
                ent_b = active[idx_b]

                if ent_a.chain != ent_b.chain:
                    continue
                if len(ent_a.crossings) != len(ent_b.crossings):
                    continue
                if not same_chirality(ent_a, ent_b):
                    continue
                if not loops_overlap(ent_a, ent_b):
                    continue

                crossing_distances = [
                    abs(crossing_number(cr_a) - crossing_number(cr_b))
                    for cr_a, cr_b in zip(ent_a.crossings, ent_b.crossings)
                ]
                if not all(distance <= 20 for distance in crossing_distances):
                    continue

                loop_a = ent_a.j - ent_a.i
                loop_b = ent_b.j - ent_b.i
                shorter = ent_a if loop_a <= loop_b else ent_b
                longer = ent_b if shorter is ent_a else ent_a
                kept = merge_into_representative(shorter, longer)

                remove_indices = sorted((idx_a, idx_b), reverse=True)
                for idx in remove_indices:
                    del active[idx]
                active.append(kept)
                changed = True
                break

    return active


def representative_for_cluster(cluster: list[Entanglement]) -> Entanglement:
    if len(cluster) == 1:
        return cluster[0]

    coords = np.asarray([entanglement_coord(ent) for ent in cluster], dtype=float)
    crossing_coords = coords[:, 2:]
    median = np.median(crossing_coords, axis=0)
    distances = np.linalg.norm(crossing_coords - median, axis=1)
    loop_lengths = np.asarray([ent.j - ent.i for ent in cluster])

    best_index = min(
        range(len(cluster)),
        key=lambda idx: (distances[idx], loop_lengths[idx], cluster[idx].i, cluster[idx].j),
    )
    rep = cluster[best_index]
    num_contacts = sum(ent.num_contacts for ent in cluster)
    contacts = tuple(dict.fromkeys(contact for ent in cluster for contact in ent.contacts))
    return replace(rep, num_contacts=num_contacts, contacts=contacts)


def spatially_cluster_entanglements(entanglements: list[Entanglement], cutoff: float) -> list[Entanglement]:
    grouped = defaultdict(list)
    for ent in entanglements:
        chirality = tuple(crossing[0] for crossing in ent.crossings)
        # Keep single-crossing, double-crossing, etc. entanglements in separate
        # clustering classes so representatives never mix crossing counts.
        grouped[(ent.chain, len(ent.crossings), chirality)].append(ent)

    representatives = []
    for group in grouped.values():
        if len(group) == 1:
            representatives.append(group[0])
            continue

        coords = np.asarray([entanglement_coord(ent) for ent in group], dtype=float)
        pairs = cKDTree(coords).query_pairs(cutoff)

        parent = list(range(len(group)))

        def find(idx):
            while parent[idx] != idx:
                parent[idx] = parent[parent[idx]]
                idx = parent[idx]
            return idx

        def union(idx_a, idx_b):
            root_a = find(idx_a)
            root_b = find(idx_b)
            if root_a != root_b:
                parent[root_b] = root_a

        for idx_a, idx_b in pairs:
            union(idx_a, idx_b)

        clusters = defaultdict(list)
        for idx, ent in enumerate(group):
            clusters[find(idx)].append(ent)

        representatives.extend(representative_for_cluster(cluster) for cluster in clusters.values())

    return representatives


def cluster_glink(csv_file: str, cutoff: float) -> pd.DataFrame:
    """
    Cluster glink.py output into representative entanglements.

    The input CSV must contain the standalone GLN columns, especially `chain`,
    `resid_i`, `resid_j`, `gn`, `gc`, `Gn`, `Gc`, and `crossings`. Clustering is
    performed separately for each crossing count, so contacts with one crossing
    cannot be merged into representatives with two crossings, and so on.
    """
    raw_entanglements = read_glink_csv(csv_file)
    if not raw_entanglements:
        return output_dataframe([])

    num_raw = len(raw_entanglements)
    entanglements = select_minimal_loop_by_crossing_set(raw_entanglements)
    entanglements = merge_larger_overlapping_loops(entanglements)
    entanglements = spatially_cluster_entanglements(entanglements, cutoff)

    num_tracked = sum(ent.num_contacts for ent in entanglements)
    if num_tracked != num_raw:
        raise ValueError(f"Tracked contacts after clustering {num_tracked} != raw entangled contacts {num_raw}")

    return output_dataframe(entanglements)


def output_dataframe(entanglements: list[Entanglement]) -> pd.DataFrame:
    rows = []
    for cluster_id, ent in enumerate(
        sorted(entanglements, key=lambda item: (item.chain, item.i, item.j, item.crossings))
    ):
        rows.append(
            {
                "cluster_id": cluster_id,
                "chain": ent.chain,
                "resid_i": ent.i,
                "resid_j": ent.j,
                "gn": round(ent.gn, 5),
                "gc": round(ent.gc, 5),
                "Gn": ent.Gn,
                "Gc": ent.Gc,
                "crossings": " ".join(ent.crossings),
                "num_contacts": ent.num_contacts,
                "contacts": ";".join(ent.contacts),
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "cluster_id",
            "chain",
            "resid_i",
            "resid_j",
            "gn",
            "gc",
            "Gn",
            "Gc",
            "crossings",
            "num_contacts",
            "contacts",
        ],
    )


def resolve_cutoff(organism: str | None, cutoff: float | None) -> float:
    if cutoff is not None:
        return cutoff
    if organism is None:
        raise ValueError("Specify either --organism or --cutoff")
    if organism not in ORGANISM_CUTOFFS:
        raise ValueError("Must specify Human, Yeast, or Ecoli")
    return ORGANISM_CUTOFFS[organism]


def default_output_path(csv_file: str, outpath: str) -> Path:
    return Path(outpath) / f"{Path(csv_file).stem}_clustered.csv"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cluster representative entanglements from glink CSV output."
    )
    parser.add_argument("-r", "--glink_csv", required=True, help="Raw entanglement CSV from glink.")
    parser.add_argument("-o", "--outpath", required=True, help="Output directory.")
    parser.add_argument("-g", "--organism", choices=sorted(ORGANISM_CUTOFFS), help="Use organism-specific clustering cutoff.")
    parser.add_argument("-c", "--cutoff", type=float, help="Override spatial clustering cutoff.")
    parser.add_argument("-w", "--output", help="Full output CSV path. Overrides --outpath default naming.")
    args = parser.parse_args()

    outpath = Path(args.outpath)
    outpath.mkdir(parents=True, exist_ok=True)
    cutoff = resolve_cutoff(args.organism, args.cutoff)
    output = Path(args.output) if args.output else default_output_path(args.glink_csv, args.outpath)

    start = time.perf_counter()
    clustered = cluster_glink(args.glink_csv, cutoff)
    output_data = clustered.copy()
    for column in ("gn", "gc"):
        output_data[column] = output_data[column].map(lambda value: f"{value:.3f}")
    output_data.to_csv(output, index=False)
    elapsed = time.perf_counter() - start

    print(f"Wrote {len(clustered)} representative entanglements to {output}")
    print(f"Runtime: {elapsed:.6f} seconds")


if __name__ == "__main__":
    main()
