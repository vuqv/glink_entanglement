# glink-entanglement

`glink-entanglement` calculates Gaussian-linking entanglement candidates from an all-atom PDB, confirms residue crossings with Topoly, and clusters similar contacts into representative entanglements.

Author: Quyen Vu

Recommended package naming:

- distribution name: `glink-entanglement`
- Python package name: `glink_entanglement`

Current command-line scripts:

- `glink`: calculate GLN values, call Topoly, and write raw entangled contacts as CSV.
- `glink-cluster`: cluster the raw CSV into representative entanglements.

Source-tree compatibility wrappers are also provided:

- `python glink.py`
- `python clustering_glink.py`

## Package Layout

```text
glink-entanglement/
├── pyproject.toml
├── README.md
├── glink.py
├── clustering_glink.py
└── glink_entanglement/
    ├── __init__.py
    ├── __main__.py
    ├── glink.py
    └── clustering.py
```

## Installation

From the repository root:

```bash
pip install .
```

For editable development:

```bash
pip install -e .
```

## Required Modules

The standalone workflow requires:

- `MDAnalysis`
- `numpy`
- `pandas`
- `scipy`
- `topoly`

## Step 1: Calculate GLN And Crossings

Script:

```bash
glink
```

Purpose:

`glink` reads an all-atom PDB file, detects heavy-atom residue contacts, calculates partial Gaussian linking values for each contact, and uses Topoly to confirm crossing residues.

### Contact Definition

A contact is defined within each chain:

- atoms are selected with `protein and not name H*`;
- two residues are a contact if any heavy-atom pair is within `4.5 A`;
- residue sequence separation must satisfy `|j - i| >= 4`.

### GLN Calculation

For each chain:

1. C-alpha atoms define the chain path.
2. Heavy atoms define contacts.
3. For each contact `(i, j)`, residues `i..j` define the loop.
4. `gn` is calculated between the loop and the N-terminal side.
5. `gc` is calculated between the loop and the C-terminal side.
6. Absolute `gn` and `gc` are rounded into `Gn` and `Gc` using the GLN threshold.

Default threshold:

```text
0.6
```

### Topoly Crossing Confirmation

GLN alone is only a candidate screen.

A contact is written to the final CSV only if:

- `Gn != 0` or `Gc != 0`; and
- Topoly finds at least one crossing residue on a side whose rounded GLN value is nonzero.

Examples:

- `Gn != 0`, `Gc == 0`: keep only if Topoly finds `crossingsN`.
- `Gn == 0`, `Gc != 0`: keep only if Topoly finds `crossingsC`.
- `Gn != 0`, `Gc != 0`: keep if Topoly finds `crossingsN` or `crossingsC`.
- no relevant Topoly crossings: drop the contact.

Topoly loop indices are passed as 1-based indices. Topoly crossing indices are mapped back to PDB residue IDs.

### Command

```bash
glink -f testset/2ww4.pdb
```

With an output directory:

```bash
glink -f testset/2ww4.pdb -o out
```

This writes:

```text
out/2ww4_glink_contacts.csv
```

With an explicit output file:

```bash
glink -f testset/2ww4.pdb -o out/OUTPUT.csv
```

With optional parameters:

```bash
glink \
  -f testset/2ww4.pdb \
  -o out/2ww4_glink_contacts.csv \
  --GLN_threshold 0.6 \
  --contact_cutoff 4.5 \
  --topoly_density 0 \
  --topoly_min_dist 10 6 5
```

### Important Topoly Note

The standalone workflow defaults to:

```text
--topoly_density 0
```

Topoly's package default is density `1`, which is slower. Use this only when you need that higher-density surface:

```bash
glink -f testset/2ww4.pdb --topoly_density 1
```

### Output

If `-o` is omitted, the default output name is:

```text
<pdb_stem>_glink_contacts.csv
```

If `-o` is a directory, the same default filename is written inside that directory. If `-o` ends in `.csv`, it is treated as the full output file path.

CSV columns:

| Column | Description |
| --- | --- |
| `chain` | Chain/segment identifier. |
| `resid_i` | PDB residue ID for the first contact residue. |
| `resname_i` | Residue name for `resid_i`. |
| `resid_j` | PDB residue ID for the second contact residue. |
| `resname_j` | Residue name for `resid_j`. |
| `contact_i_index` | Zero-based C-alpha index for the first contact residue. |
| `contact_j_index` | Zero-based C-alpha index for the second contact residue. |
| `gn` | Raw N-terminal partial Gaussian linking value. |
| `gc` | Raw C-terminal partial Gaussian linking value. |
| `Gn` | Rounded absolute `gn`. |
| `Gc` | Rounded absolute `gc`. |
| `crossingsN` | Topoly N-terminal crossing residues mapped to PDB residue IDs. |
| `crossingsC` | Topoly C-terminal crossing residues mapped to PDB residue IDs. |
| `crossings` | Combined crossing residues used by clustering. |

The script prints total runtime:

```text
Runtime: <seconds> seconds
```

## Step 2: Cluster Representative Entanglements

Script:

```bash
glink-cluster
```

Purpose:

`glink-cluster` reads the CSV from `glink` and clusters similar entangled contacts into representative entanglements.

### Clustering Logic

The script performs these stages:

1. Read raw entangled contacts from the `glink` CSV.
2. Drop rows without crossing residues.
3. Group by chain and exact crossing set.
4. For each exact crossing set, keep the shortest loop as the first representative and track all raw contacts represented.
5. Merge larger overlapping loops only when they have:
   - the same chain;
   - the same number of crossings;
   - the same chirality sequence;
   - overlapping loop endpoints;
   - corresponding crossing residues within 20 residues.
6. Spatially cluster representatives by Euclidean distance over:

```text
(resid_i, resid_j, crossing_residue_ids...)
```

Clustering is separated by:

```text
(chain, number_of_crossings, chirality_sequence)
```

This means one-crossing representatives are not merged into two-crossing representatives, and so on.

### Organism Cutoffs

Preset cutoffs:

| Organism | Cutoff |
| --- | --- |
| `Human` | `52` |
| `Ecoli` | `57` |
| `Yeast` | `49` |

You can also provide a custom cutoff with `--cutoff`.

### Command

Using an organism preset:

```bash
glink-cluster \
  -r 2ww4_glink_contacts.csv \
  -o clustered_glink \
  -g Human
```

Using a custom cutoff:

```bash
glink-cluster \
  -r 2ww4_glink_contacts.csv \
  -o clustered_glink \
  -c 52
```

Writing to an explicit output file:

```bash
glink-cluster \
  -r 2ww4_glink_contacts.csv \
  -o clustered_glink \
  -g Human \
  -w clustered_glink/2ww4_representative_entanglements.csv
```

### Arguments

| Short | Long | Description |
| --- | --- | --- |
| `-r` | `--glink_csv` | Raw entanglement CSV from `glink`. |
| `-o` | `--outpath` | Output directory. |
| `-g` | `--organism` | Organism preset: `Human`, `Ecoli`, or `Yeast`. |
| `-c` | `--cutoff` | Custom spatial clustering cutoff. |
| `-w` | `--output` | Full output CSV path. |

### Output

Default output name:

```text
<input_csv_stem>_clustered.csv
```

CSV columns:

| Column | Description |
| --- | --- |
| `cluster_id` | Integer representative cluster ID. |
| `chain` | Chain/segment identifier. |
| `resid_i` | Representative loop start residue ID. |
| `resid_j` | Representative loop end residue ID. |
| `gn` | Raw `gn` value from the representative contact. |
| `gc` | Raw `gc` value from the representative contact. |
| `Gn` | Rounded absolute `gn` from the representative contact. |
| `Gc` | Rounded absolute `gc` from the representative contact. |
| `crossings` | Representative crossing residue set. |
| `num_contacts` | Number of raw contacts represented by the cluster. |
| `contacts` | Semicolon-delimited raw contact list, formatted as `i-j`. |

## End-To-End Example

```bash
glink \
  -f testset/2ww4.pdb \
  -o out

glink-cluster \
  -r out/2ww4_glink_contacts.csv \
  -o clustered_glink \
  -g Human
```

## Notes

- `glink_entanglement.glink` uses `pathlib.Path` for output path handling.
- `glink_entanglement.clustering` also uses `pathlib.Path` and creates the output directory if needed.
- The workflow is CSV-based.
