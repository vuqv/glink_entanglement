# Gaussian Entanglement Pipeline Specifications

This repository contains a two-stage pipeline for identifying and reducing protein Gaussian linking number entanglements.

1. `gaussian_entanglement.py` computes raw, unmapped entanglements from PDB structures.
2. `clustering.py` clusters the raw entanglements into representative entanglements.

The specifications below describe the behavior implemented in the current code.

## Stage 1: Raw Gaussian Entanglement Calculation

### Script

`gaussian_entanglement.py`

### Purpose

Calculate raw non-covalent lasso-like entanglements for one PDB file or all PDB files in a directory. Each entanglement is defined by:

- a native-contact loop, represented as residue pair `(i, j)`;
- one or more crossing residues, each labeled with chirality, such as `+123` or `-123`;
- partial Gaussian linking values against the N-terminal and C-terminal chain regions, `gn` and `gc`;
- a disulfide-bond flag for the native-contact residue pair.

### Command Line Interface

```bash
python gaussian_entanglement.py \
  --PDB <pdb-file-or-directory> \
  [--GLN_threshold <float>] \
  [--Calpha True|False] \
  [--topoly_density <int>]
```

### Arguments

| Argument | Required | Default | Description |
| --- | --- | --- | --- |
| `--PDB` | Yes | None | Path to a PDB file or a directory containing PDB files. |
| `--GLN_threshold` | No | `0.6` | Fractional threshold used by `point_rounding()` to round absolute `gn`/`gc` values to an integer entanglement count. |
| `--Calpha` | No | `False` | If set to string `True`, use C-alpha atoms and an 8 Angstrom contact cutoff. Otherwise use non-hydrogen atoms and a 4.5 Angstrom contact cutoff. |
| `--topoly_density` | No | `0` | Density passed to `topoly.lasso_type()` for minimal-loop surface triangulation. Higher values may reduce unrealistic crossings for disordered loops at higher runtime cost. |

### Dependencies

The script imports:

- `Bio.PDB.PDBParser`
- `MDAnalysis`
- `numpy`
- `pandas`
- `numba`
- `scipy.spatial.distance`
- `topoly.lasso_type`

`Bio.PDB.DSSP`, `multiprocessing`, and `sys` are imported but are not used by the active calculation path.

### Created Directories

On startup, the script creates missing output directories in the current working directory:

| Mode | Directories |
| --- | --- |
| default heavy-atom mode | `unmapped_GE`, `unmapped_missing_residues`, `Before_TER_PDBs`, `clustered_unmapped_GE` |
| C-alpha mode | `unmapped_GE_CA`, `unmapped_missing_residues`, `Before_TER_PDBs`, `clustered_unmapped_GE` |

`Before_TER_PDBs` is created, but the current CLI does not call `pre_processing_pdb()`.

### Input Handling

If `--PDB` is a file, that file is processed.

If `--PDB` is a directory, the script builds input paths as:

```python
f"{pdb_dir}{filename}"
```

Therefore, directory inputs are expected to include a trailing path separator, for example `PDBs/`.

Each PDB is loaded with:

```python
mda.Universe(pdb_file, format="PDB")
```

The protein identifier used in output names is the input file basename without its extension.

### Per-PDB Workflow

For each PDB:

1. Skip the PDB if its raw output file already exists in the relevant raw output directory.
2. Load the PDB with MDAnalysis.
3. Detect disulfide bonds by scanning cysteine `SG` atom pairs with distance `< 2.2`.
4. Select atoms:
   - default mode: `not name H*`;
   - C-alpha mode: `name CA`.
5. Iterate over MDAnalysis segment IDs as chains.
6. For each chain:
   - remove duplicate `(resid, atom name)` rows;
   - build residue-index to atom-index mappings;
   - require at least two C-alpha atoms;
   - verify mapped coordinates against MDAnalysis selections;
   - compute raw entanglements;
   - write raw entanglements and missing residues.

### Native Contact Definition

Native contacts are built from pairwise distances:

| Mode | Atom selection | Contact cutoff |
| --- | --- | --- |
| default | non-hydrogen atoms | `4.5` Angstrom |
| C-alpha | C-alpha atoms | `8.0` Angstrom |

Contacts are restricted to the upper triangle with offset `k=4`, then reduced to residue-level contacts. A residue pair is retained only if the absolute PDB residue-number separation is greater than 4.

### Gaussian Linking Calculation

For a chain with `l` residues:

1. C-alpha midpoints and segment vectors are built for consecutive residue pairs.
2. A full `(l - 1) x (l - 1)` dot-product matrix is computed from the Gauss linking integral terms.
3. For each native-contact loop `(i, j)`:
   - loop range is `i` through `j - 1`;
   - N-terminal comparison range is `5` through `i - 6`;
   - C-terminal comparison range is `j + 6` through `l - 6`;
   - `gn` and `gc` are sums of the relevant matrix entries divided by `4*pi`.
4. Absolute `gn` and `gc` values are rounded by `point_rounding()` using `GLN_threshold`.
5. A native contact is considered a candidate entanglement if either rounded value has absolute value at least 1.

### Crossing Detection

Crossing residues are determined with:

```python
topoly.lasso_type(
    coor,
    loop_indices=native_contacts,
    more_info=True,
    density=density,
    min_dist=[10, 6, 5],
)
```

The implemented comments define the reduction criteria as:

- crossing residues must be at least 10 residues apart;
- the first crossing must be at least 6 residues from the loop;
- the first crossing must be at least 5 residues from the nearest terminus.

For each candidate loop:

- N-terminal crossings are included if rounded `gn >= 1` in absolute value.
- C-terminal crossings are included if rounded `gc >= 1` in absolute value.
- Crossing labels are converted from Topoly residue indices to PDB residue numbers while preserving chirality.

### Missing-Residue Filters

Missing residues are inferred as gaps between the first and last C-alpha residue IDs.

Loop candidates are removed if either:

- the loop contains at least three consecutive missing residues; or
- missing residues exceed 5% of the loop length.

Final entanglements are removed if any missing residue lies within plus/minus 10 residues of any crossing residue.

### Raw Output Files

Raw entanglement output is appended to:

| Mode | Output path |
| --- | --- |
| default | `unmapped_GE/<protein>_GE.txt` |
| C-alpha | `unmapped_GE_CA/<protein>_GE.txt` |

Each raw entanglement line has this pipe-delimited schema:

```text
Chain <chain> | (<i>, <j>, <crossings_array>) | <gn> | <gc> | CCbond-<True|False>
```

Example shape:

```text
Chain A | (10, 85, ['+42' '-57']) | 1.23456 | -0.12345 | CCbond-False
```

Notes:

- `<i>` and `<j>` are PDB residue IDs, not zero-based residue indices.
- `<crossings_array>` is written from a NumPy unique array, so entries are space-separated inside brackets.
- `gn` and `gc` are raw partial Gaussian linking values, not rounded values.
- `CCbond-True` means the native-contact pair was detected as a cysteine disulfide bond.

### Missing Residue Output

If missing residues are found, they are appended to:

```text
unmapped_missing_residues/<protein>_M.txt
```

Each line has this format:

```text
Chain <chain>: <resid_1> <resid_2> ...
```

### Termination

After all inputs are processed, the script prints:

```text
NORMAL TERMINATION
```

## Stage 2: Representative Entanglement Clustering

### Script

`clustering.py`

### Purpose

Read one raw entanglement file from `gaussian_entanglement.py`, merge and spatially cluster related raw entanglements, then write representative entanglements with the number and list of raw native contacts represented by each cluster.

### Command Line Interface

```bash
python clustering.py \
  --prot_unmapped_GE_file <raw-ge-file> \
  --outpath <output-directory> \
  --organism Human|Ecoli|Yeast
```

### Arguments

| Argument | Required | Description |
| --- | --- | --- |
| `--prot_unmapped_GE_file` | Yes | Path to a raw entanglement file produced by `gaussian_entanglement.py`. |
| `-o`, `--outpath` | Yes | Directory where clustered output should be written. Created if missing. |
| `--organism` | Yes | Selects the spatial clustering distance cutoff. Must be `Human`, `Ecoli`, or `Yeast`. |

### Organism-Specific Cutoffs

| Organism | Distance cutoff |
| --- | --- |
| `Human` | `52` |
| `Ecoli` | `57` |
| `Yeast` | `49` |

Any other organism value raises `ValueError("Must specify Human, Yeast, or Ecoli")`.

### Dependencies

The script imports:

- `numpy`
- `geom_median.numpy.compute_geometric_median`
- `scipy.spatial.distance.cdist`

It also uses standard-library modules `collections`, `itertools`, `functools`, `re`, `math`, `random`, `copy`, and `os`.

### Input File Requirements

The input file must use the raw output format from `gaussian_entanglement.py`.

Only lines that split into exactly five pipe-delimited fields are processed:

```text
chain | ijr | gn | gc | CCbond
```

The current parser expects:

- the filename to contain the protein ID before the first underscore;
- raw filenames to contain only one underscore in normal use, otherwise a warning is printed;
- the `ijr` field to contain exactly two native-contact residues followed by one crossing-residue array;
- crossing residues to include a leading chirality sign, such as `+42` or `-42`.

### Output File

The clustered output file is:

```text
<outpath>/<protein>_clustered_GE.txt
```

If the file already exists, clustering is skipped and the script asks the user to delete the file before reclustering.

The output header is:

```text
chain|ijr|gn|gc|num_contacts|contacts|CCBond
```

Each output row has this schema:

```text
<chain>|(<i>, <j>, <crossing_1>, ...)|<gn>|<gc>|<num_contacts>|<contacts>|<CCBond>
```

Where:

- `chain` is the chain label from the raw file, for example `Chain A`.
- `ijr` is the representative native-contact pair plus crossing residues.
- `gn` and `gc` are copied from the representative raw entanglement and formatted to five decimals.
- `num_contacts` is the number of raw native contacts represented by the output row.
- `contacts` is a semicolon-delimited list of raw contact pairs, each formatted as `i-j`.
- `CCBond` is `True` if any represented contact is a detected disulfide bond, otherwise `False`.

### Clustering Workflow

`cluster_entanglements()` performs five implemented stages.

#### Step 1: Parse and Group Raw Entanglements

For each raw line:

1. Parse chain, native-contact pair, crossings, `gn`, `gc`, and `CCbond`.
2. Sort crossing residues by numeric residue ID while preserving chirality.
3. Group native-contact loops by `(chain, sorted crossings)`.
4. Store `gn` and `gc` by `(i, j, crossings...)`.
5. Track the starting number of raw entanglements per chain for QC.
6. Track raw contacts marked `CCbond-True`.

#### Step 2a: Select Minimal Loop Per Crossing Set

For each `(chain, crossing set)` group:

1. Compute the length of each loop as `j - i`.
2. Select the shortest loop as the representative loop for that crossing set.
3. Store an entanglement tuple:

```text
(num_loops, representative_i, representative_j, crossings..., contacts)
```

`contacts` is a semicolon-delimited list of all raw loops in the group.

The code checks that the summed `num_loops` still equals the raw entanglement count for each chain.

#### Step 2b: Merge Nearby Different-Size Crossing Sets

Pairs of entanglements are considered for merging when:

- at least one crossing pair has the same chirality and residue distance `<= 3`;
- the two entanglements have different numbers of crossing residues;
- either endpoint of one loop lies within the inclusive range of the other loop;
- no crossing residue lies inside the combined loop range from `min(i, j, k, l)` to `max(i, j, k, l)`;
- the minimum assignment distance between the smaller crossing set and the larger crossing set is `<= 20`.

When these criteria pass, the entanglement with more crossing residues is kept and the smaller one is merged into it. The raw-contact count and contact list are updated so the total number of raw entanglements remains conserved.

If merged records cannot be assigned back to a kept representative, the code raises a `ValueError` that suggests deleting the problematic entanglement or increasing Topoly density in the raw calculation.

#### Step 3: Remove Larger Overlapping Loops

Pairs of processed entanglements are considered when:

- they have the same number of crossing residues;
- their chirality sequences are identical;
- either endpoint of one loop lies within the inclusive range of the other loop;
- corresponding crossing residues are all within `20` residues of each other.

When these criteria pass, the larger loop is removed and merged into the smaller loop. If loop lengths are equal, the first entanglement is treated as the larger loop and the second as the shorter loop.

The code again checks that the summed raw-contact counts remain equal to the starting raw count per chain.

#### Step 4: Spatial Clustering

Processed entanglements are grouped by:

```text
<chain>_<number_of_crossings>_<chirality_sequence>
```

Within each group:

1. Pairwise distances are computed by `loop_distance()`.
2. `loop_distance()` removes crossing chiralities and computes Euclidean distance over:

```text
(i, j, crossing_residue_ids...)
```

3. Pairs with distance `<= organism_cutoff` are linked into clusters, with duplicate captures suppressed.
4. Unclustered entanglements become singleton clusters.
5. For multi-entanglement clusters:
   - compute the geometric median of crossing residue positions;
   - choose candidates with minimum distance to that median;
   - choose the candidate with the shortest loop length;
   - break remaining ties randomly;
   - aggregate `num_contacts` and semicolon-delimited contact lists.
6. Singleton clusters are retained directly.

The final QC check verifies that the sum of `num_contacts` in clustered representatives equals the raw entanglement count for every chain.

#### Step 5: Write Clustered Representatives

The output file is written with one representative row per final cluster.

For each representative:

- `gn` and `gc` are looked up from the raw representative `(i, j, crossings...)`;
- `num_contacts` and `contacts` represent all raw contacts merged into the cluster;
- `CCBond` is set by checking whether any stored disulfide contact appears in the represented contact list.

### Determinism

The clustering stage uses `random.choice()` in two places:

- when selecting among cluster seeds with the same number of neighbors;
- when selecting among representative candidates tied by geometric-median distance and loop length.

No random seed is set. Therefore, tied clustering results may not be deterministic across runs.

### Quality Checks and Failure Modes

The clustering code raises errors if raw entanglement counts are not conserved after:

- Step 2a;
- Step 2b;
- Step 3;
- final representative clustering.

It may also raise parsing or lookup errors when the raw file does not match the expected format or when representative `gn`/`gc` data cannot be found.

## End-to-End Usage

Default heavy-atom workflow:

```bash
python gaussian_entanglement.py --PDB PDBs/example.pdb
python clustering.py \
  --prot_unmapped_GE_file unmapped_GE/example_GE.txt \
  --outpath clustered_unmapped_GE \
  --organism Human
```

C-alpha workflow:

```bash
python gaussian_entanglement.py --PDB PDBs/example.pdb --Calpha True
python clustering.py \
  --prot_unmapped_GE_file unmapped_GE_CA/example_GE.txt \
  --outpath clustered_unmapped_GE \
  --organism Human
```

## Important Implementation Assumptions

- Chain IDs come from MDAnalysis segment IDs.
- Residue IDs are treated as the biological residue labels from the PDB, not normalized zero-based sequence positions.
- Directory input to `gaussian_entanglement.py` should include a trailing slash.
- Existing output files cause the corresponding calculation or clustering step to be skipped.
- Raw output is append-based, so deleting stale output files is required before recomputing a protein from scratch.
- `clustering.py` expects raw data produced by this exact implementation; it is not a general parser for arbitrary entanglement tables.
- `clustering.py` derives the protein ID from the first underscore-separated token in the raw filename.
- `os.sched_getaffinity(0)` is called in both scripts; environments that do not provide this function may fail before calculation starts.

