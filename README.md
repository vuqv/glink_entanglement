# glink-entanglement

`glink-entanglement` calculates Gaussian-linking entanglement candidates from an all-atom PDB, confirms residue crossings with Topoly, and clusters similar contacts into representative entanglements.

Author: Quyen Vu

Recommended package naming:

- distribution name: `glink-entanglement`
- Python package name: `glink_entanglement`

Current command-line scripts:

- `glink`: calculate GLN values, call Topoly, and write raw entangled contacts as CSV.
- `glink-cluster`: cluster the raw CSV into representative entanglements while treating `crossingsN` and `crossingsC` as different crossing sets.

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

Topoly accepts 1-based coordinate-array indices. `glink` stores contact indices as zero-based C-alpha indices in the CSV columns `i` and `j`, adds `1` before passing loop indices to Topoly, and treats Topoly crossing labels as 1-based indices when mapping them back to PDB residue IDs.

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
| `contact` | Contact label formatted as `<resname_i><resid_i>-<resname_j><resid_j>`, for example `ILE17-THR194`. |
| `i` | Zero-based C-alpha index for the first contact residue. |
| `j` | Zero-based C-alpha index for the second contact residue. |
| `gn` | N-terminal partial Gaussian linking value, printed to three decimals. |
| `gc` | C-terminal partial Gaussian linking value, printed to three decimals. |
| `Gn` | Rounded absolute `gn`. |
| `Gc` | Rounded absolute `gc`. |
| `crossingsN` | Topoly N-terminal crossing residues. Topoly reports 1-based crossing indices; the CSV values are mapped to PDB residue IDs. |
| `crossingsC` | Topoly C-terminal crossing residues. Topoly reports 1-based crossing indices; the CSV values are mapped to PDB residue IDs. |
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

Clustering keeps N-terminal and C-terminal Topoly crossings separate. For example, `crossingsN = +45` and `crossingsC = +45` are treated as different crossing signatures and cannot be merged into the same representative just because the residue number is the same.

Clustering is separated by:

```text
(chain, number_of_crossings, side_and_chirality_sequence)
```

This means one-crossing representatives are not merged into two-crossing representatives, and N-side crossings are not merged with C-side crossings.

### Organism Cutoffs

Preset cutoffs:

| Organism | Cutoff |
| --- | --- |
| `Custom` | `52` |
| `Human` | `52` |
| `Ecoli` | `57` |
| `Yeast` | `49` |

`Custom` is the default when `-g/--organism` is not specified. You can also provide a custom cutoff with `--cutoff`; this overrides the organism preset.

### Command

Using the default `Custom` cutoff of 52:

```bash
glink-cluster \
  -f 2ww4_glink_contacts.csv \
  -o clustered_glink
```

Using an organism preset:

```bash
glink-cluster \
  -f 2ww4_glink_contacts.csv \
  -o clustered_glink \
  -g Human
```

Using an explicit custom cutoff:

```bash
glink-cluster \
  -f 2ww4_glink_contacts.csv \
  -o clustered_glink \
  -c 52
```

Writing to an explicit output file:

```bash
glink-cluster \
  -f 2ww4_glink_contacts.csv \
  -o clustered_glink/2ww4_representative_entanglements.csv \
  -g Human
```

### Arguments

| Short | Long | Description |
| --- | --- | --- |
| `-f` | `--glink_csv` | Required. Raw entanglement CSV from `glink`. |
| `-o` | `--output_path`, `--outpath` | Required. Output directory or CSV path. If a directory is provided, writes `<input_csv_stem>_clustered.csv` inside it. If a `.csv` path is provided, writes exactly to that file. |
| `-g` | `--organism` | Optional organism preset: `Custom`, `Human`, `Ecoli`, or `Yeast`. Default: `Custom`. |
| `-c` | `--cutoff` | Optional custom spatial clustering cutoff. Overrides organism preset. |
| `-w` | `--output` | Optional legacy explicit output CSV path. Overrides `-o/--output_path`. |

### Output

Default output behavior:

```text
-o clustered_glink      -> clustered_glink/<input_csv_stem>_clustered.csv
-o clustered/file.csv   -> clustered/file.csv
```

CSV columns:

| Column | Description |
| --- | --- |
| `cluster_id` | Integer representative cluster ID. |
| `contact` | Representative contact label formatted as `<resname_i><resid_i>-<resname_j><resid_j>`. |
| `i` | Zero-based C-alpha index for the representative contact's first residue. |
| `j` | Zero-based C-alpha index for the representative contact's second residue. |
| `gn` | `gn` value from the representative contact, printed to three decimals. |
| `gc` | `gc` value from the representative contact, printed to three decimals. |
| `Gn` | Rounded absolute `gn` from the representative contact. |
| `Gc` | Rounded absolute `gc` from the representative contact. |
| `crossingsN` | Representative N-terminal crossing residue set. |
| `crossingsC` | Representative C-terminal crossing residue set. |
| `crossings` | Combined `crossingsN` and `crossingsC`, retained for readability and compatibility. |
| `num_contacts` | Number of raw contacts represented by the cluster. |
| `contacts` | Semicolon-delimited raw contact list, formatted like `ILE17-THR194` when the input CSV uses contact labels. |

## End-To-End Example

```bash
glink \
  -f testset/2ww4.pdb \
  -o out

glink-cluster \
  -f out/2ww4_glink_contacts.csv \
  -o clustered_glink \
  -g Human
```

## Notes

- `glink_entanglement.glink` uses `pathlib.Path` for output path handling.
- `glink_entanglement.clustering` also uses `pathlib.Path` and creates the output directory if needed.
- The workflow is CSV-based.
