---
layout: post
title: GLink Entanglement Tutorial
date: 2026-05-31 12:00:00-0400
description: A practical tutorial for detecting, confirming, and clustering Gaussian-linking entanglements in protein structures.
tags: protein-entanglement gaussian-linking topology python structural-biology
categories: research-notes
giscus_comments: true
related_posts: false
---
`glink-entanglement` is a Python workflow for finding candidate loop-threading entanglements in all-atom protein structures. It starts from heavy-atom residue contacts, scores each contact with partial Gaussian linking values, confirms the retained candidates with Topoly crossing residues, and then clusters many contact-level hits into representative entanglements.

The package currently exposes two command-line programs:

- `glink` calculates contact-level Gaussian-linking entanglements from a PDB file.
- `glink-cluster` converts the raw contact table into a smaller representative table.

The workflow is CSV based. The raw table is the audit trail; the clustered table is usually the table to compare across structures.

## When to use it

Use this workflow when you have an all-atom protein PDB and want a residue-level screen for native-contact loop entanglements. Typical inputs include experimental PDB structures, AlphaFold-style PDB models, or representative structures extracted from a molecular dynamics trajectory.

The method analyzes protein chains separately. For each chain, the package selects:

```text
protein and not name H*
```

Hydrogens are ignored for contact detection. C-alpha atoms define the curve used for the Gaussian-linking calculation and for Topoly.

## Installation

Clone the package and install it in the active Python environment:

```bash
git clone https://github.com/vuqv/glink_entanglement.git
cd glink_entanglement
pip install .
```

For development, install in editable mode:

```bash
pip install -e .
```

The package requires Python 3.10 or newer and depends on `MDAnalysis`, `numpy`, `pandas`, `scipy`, and `topoly`.

After installation, confirm that both scripts are on your path:

```bash
glink --help
glink-cluster --help
```

The source tree also includes compatibility wrappers:

```bash
python glink.py --help
python clustering_glink.py --help
```

## The two-stage workflow

A typical analysis has two commands:

```bash
mkdir -p out clustered

glink \
  -f PATH_TO_PDB \
  -o out

glink-cluster \
  -f out/PDB_STEM_glink_contacts.csv \
  -o clustered
```

The first command writes a raw contact CSV named from the input PDB stem:

```text
out/PDB_STEM_glink_contacts.csv
```

The second command writes a clustered representative CSV:

```text
clustered/PDB_STEM_glink_contacts_clustered.csv
```

If `-o` points to a directory, the package creates that directory and writes the default filename inside it. If `-o` ends in `.csv`, the package treats it as the exact output file.

## Stage 1: raw entangled contacts

Run `glink` on one PDB file:

```bash
glink -f PATH_TO_PDB -o out
```

The full set of tunable parameters is:

```bash
glink \
  -f PATH_TO_PDB \
  -o out/PDB_STEM_glink_contacts.csv \
  --GLN_threshold 0.6 \
  --contact_cutoff 4.5 \
  --topoly_density 0 \
  --topoly_min_dist 10 6 5
```

The defaults are a good first pass:

| Option | Default | Meaning |
| --- | --- | --- |
| `-f`, `--PDB` | required | Input all-atom PDB file. |
| `-o`, `--output` | `<pdb_stem>_glink_contacts.csv` | Output directory or explicit CSV path. |
| `--contact_cutoff` | `4.5` | Heavy-atom distance cutoff, in Angstrom, for residue contacts. |
| `--GLN_threshold` | `0.6` | Threshold for rounding absolute `gn` and `gc` into `Gn` and `Gc`. |
| `--topoly_density` | `0` | Topoly minimal-surface triangulation density. |
| `--topoly_min_dist` | `10 6 5` | Topoly crossing-reduction distances for crossing, loop, and tail end. |

At the end of a successful run, `glink` reports how many contacts were written and the runtime:

```text
Wrote N contacts to out/PDB_STEM_glink_contacts.csv
Runtime: 12.345678 seconds
```

The exact count and runtime depend on protein size, contact density, and Topoly settings.

## How contacts become candidates

Within each chain, two residues form a candidate contact when any pair of their heavy atoms is within the contact cutoff and the residues are separated by at least four positions along the sequence.

For a contact `(i, j)`, the residues from `i` through `j` define the loop. The package calculates two partial Gaussian-linking values:

| Value | Meaning |
| --- | --- |
| `gn` | Gaussian linking between the loop and the N-terminal side of the chain. |
| `gc` | Gaussian linking between the loop and the C-terminal side of the chain. |
| `Gn` | Rounded absolute value of `gn`. |
| `Gc` | Rounded absolute value of `gc`. |

Contacts with `Gn == 0` and `Gc == 0` are discarded before Topoly confirmation.

## Topoly confirmation

Gaussian linking is used as a fast screen, but it is not the final retention criterion. A contact is written only when Topoly finds at least one crossing on a side whose rounded Gaussian-linking score is nonzero:

| Rounded GLN state | Retention rule |
| --- | --- |
| `Gn != 0`, `Gc == 0` | Keep only if Topoly reports `crossingsN`. |
| `Gn == 0`, `Gc != 0` | Keep only if Topoly reports `crossingsC`. |
| `Gn != 0`, `Gc != 0` | Keep if either `crossingsN` or `crossingsC` is present. |
| no relevant crossing | Drop the contact. |

Topoly uses 1-based coordinate-array indices for loop and crossing labels. `glink` stores the public contact indices as zero-based C-alpha indices in the CSV columns `i` and `j`, adds one before passing loop endpoints to Topoly, and maps Topoly crossing labels back to PDB residue IDs for `crossingsN` and `crossingsC`.

This distinction matters when debugging: `i` and `j` are zero-based C-alpha positions in the analyzed chain, while crossing labels such as `-7` or `+172` are signed PDB residue IDs after Topoly's crossing labels have been mapped back to the input residues.

The package default is:

```text
--topoly_density 0
```

Topoly's own package default is density `1`, which is slower. Use density `1` when the denser surface calculation is worth the extra cost:

```bash
glink -f PATH_TO_PDB -o out --topoly_density 1
```

## Raw CSV output

The current command-line raw CSV has one row per retained entangled contact and writes only the public contact-level fields:

```text
contact,i,j,gn,gc,Gn,Gc,crossingsN,crossingsC,crossings
ALA10-GLY42,9,41,-0.809,-0.571,1,0,-7,,-7
SER15-THR80,14,79,0.120,-0.923,0,1,,+95,+95
```

The columns are:

| Column | Meaning |
| --- | --- |
| `contact` | Loop-closing contact label formatted as `<resname_i><resid_i>-<resname_j><resid_j>`, for example `ALA10-GLY42`. |
| `i` | Zero-based C-alpha index for the first contact residue. |
| `j` | Zero-based C-alpha index for the second contact residue. |
| `gn` | N-terminal partial Gaussian-linking value, printed to three decimals by the CLI. |
| `gc` | C-terminal partial Gaussian-linking value, printed to three decimals by the CLI. |
| `Gn` | Rounded absolute `gn`. |
| `Gc` | Rounded absolute `gc`. |
| `crossingsN` | Signed Topoly N-terminal crossing residues, mapped to PDB residue IDs. |
| `crossingsC` | Signed Topoly C-terminal crossing residues, mapped to PDB residue IDs. |
| `crossings` | Combined crossing list retained for readability and clustering compatibility. |

The raw CLI output does not include separate `chain`, `resid_i`, `resname_i`, `resid_j`, `resname_j`, `contact_i_index`, or `contact_j_index` columns. Those fields exist inside the Python-level DataFrame, but the command-line output collapses them into the readable `contact` label and the zero-based `i` and `j` indices.

The sign on a crossing is Topoly's chirality. The number is the mapped PDB residue ID, not the zero-based `i` or `j` index.

## Stage 2: clustering representative entanglements

Raw contacts can be redundant because many nearby loop-closing contacts may describe the same entanglement. `glink-cluster` groups those raw rows and reports representative entanglements:

```bash
glink-cluster \
  -f out/PDB_STEM_glink_contacts.csv \
  -o clustered
```

The current clustering input flag is `-f` or `--glink_csv`. If `-g/--organism` and `-c/--cutoff` are both omitted, the command uses the `Custom` preset, which has cutoff `52`.

| Option | Default | Meaning |
| --- | --- | --- |
| `-f`, `--glink_csv` | required | Raw CSV produced by `glink`. |
| `-o`, `--output_path`, `--outpath` | required | Output directory or explicit CSV path. |
| `-g`, `--organism` | `Custom` | Organism cutoff preset: `Custom`, `Human`, `Ecoli`, or `Yeast`. |
| `-c`, `--cutoff` | none | Explicit spatial clustering cutoff. Overrides the organism preset. |
| `-w`, `--output` | none | Legacy explicit output CSV path. Overrides `-o`. |

Preset cutoffs are:

| Organism | Cutoff |
| --- | --- |
| `Custom` | `52` |
| `Human` | `52` |
| `Ecoli` | `57` |
| `Yeast` | `49` |

The default command is therefore equivalent to:

```bash
glink-cluster \
  -f out/PDB_STEM_glink_contacts.csv \
  -o clustered \
  -g Custom
```

Use `--cutoff` when a preset is not appropriate:

```bash
glink-cluster \
  -f out/PDB_STEM_glink_contacts.csv \
  -o clustered \
  --cutoff 52
```

To write to an explicit file, pass a `.csv` path to `-o`:

```bash
glink-cluster \
  -f out/PDB_STEM_glink_contacts.csv \
  -o clustered/PDB_STEM_representative_entanglements.csv
```

The legacy `-w` option is still accepted for explicit output paths:

```bash
glink-cluster \
  -f out/PDB_STEM_glink_contacts.csv \
  -o clustered \
  -w clustered/PDB_STEM_representative_entanglements.csv
```

## How clustering works

The clustering code keeps N-terminal and C-terminal crossings separate. A crossing at residue `+45` on the N-terminal side and a crossing at residue `+45` on the C-terminal side are different signatures.

The clustering stages are:

1. Read the raw `glink` CSV.
2. Drop rows without `crossingsN` or `crossingsC`.
3. Parse the residue IDs from `contact`.
4. Group rows by chain and exact side-specific crossing set.
5. Within each exact crossing set, choose the shortest loop as the first representative.
6. Merge larger overlapping loops when they share chain, crossing count, side-specific chirality sequence, overlapping loop endpoints, and crossing residues within 20 residues.
7. Spatially cluster remaining representatives by Euclidean distance over `(residue_i, residue_j, crossing_residue_ids...)`.
8. Choose a representative near the median crossing coordinate, preferring shorter loops when tied.

The clustering output also tracks how many raw contacts were represented. If the number of tracked contacts does not match the number of raw entangled contacts, the code raises an error instead of silently dropping rows.

## Clustered CSV schema

The current source writes one row per representative entanglement:

```text
cluster_id,contact,i,j,gn,gc,Gn,Gc,crossingsN,crossingsC,crossings,num_contacts,contacts
0,ALA10-GLY42,9,41,-0.809,-0.571,1,0,-7,,-7,3,ALA10-GLY42;ALA10-SER43;VAL11-GLY42
```

The columns are:

| Column | Meaning |
| --- | --- |
| `cluster_id` | Integer representative cluster ID. |
| `contact` | Representative loop-closing contact label. |
| `i` | Zero-based C-alpha index for the representative contact's first residue. |
| `j` | Zero-based C-alpha index for the representative contact's second residue. |
| `gn` | `gn` value from the representative contact, printed to three decimals by the CLI. |
| `gc` | `gc` value from the representative contact, printed to three decimals by the CLI. |
| `Gn` | Rounded absolute `gn` from the representative contact. |
| `Gc` | Rounded absolute `gc` from the representative contact. |
| `crossingsN` | Representative N-terminal crossing set. |
| `crossingsC` | Representative C-terminal crossing set. |
| `crossings` | Combined crossing list for readability and compatibility. |
| `num_contacts` | Number of raw contacts represented by the cluster. |
| `contacts` | Semicolon-delimited raw contact labels represented by the cluster. |

For reporting, `contact`, `crossingsN`, `crossingsC`, and `num_contacts` are usually the most useful columns. For auditing, use `contacts` to return to the raw CSV rows that contributed to a representative.

## Optional: remove slipknots

After clustering, you can optionally remove slipknot crossing pairs from the clustered CSV. This is a post-processing step, not a required part of the core `glink` and `glink-cluster` workflow.

The single-file script is:

```bash
python scripts/remove_slipknots.py \
  -i clustered/PDB_STEM_glink_contacts_clustered.csv
```

By default, it writes a sibling file:

```text
clustered/PDB_STEM_glink_contacts_clustered_no_slipknot.csv
```

You can also choose an output directory or exact output file:

```bash
python scripts/remove_slipknots.py \
  -i clustered/PDB_STEM_glink_contacts_clustered.csv \
  -o no_slipknot
```

```bash
python scripts/remove_slipknots.py \
  -i clustered/PDB_STEM_glink_contacts_clustered.csv \
  -o no_slipknot/PDB_STEM_no_slipknot.csv
```

The script reads `crossingsN`, `crossingsC`, and `crossings`, then removes adjacent opposite-sign crossing pairs independently for the N-terminal and C-terminal sides. N-terminal crossings are ordered from larger residue number to smaller residue number; C-terminal crossings are ordered from smaller residue number to larger residue number. Rows whose crossings all cancel are dropped by default.

To keep rows even when all crossings cancel, use:

```bash
python scripts/remove_slipknots.py \
  -i clustered/PDB_STEM_glink_contacts_clustered.csv \
  --keep_empty
```

For many clustered CSVs, use the serial wrapper:

```bash
python scripts/run_remove_slipknots.py \
  -i batch_results/clustered \
  -o batch_results/no_slipknot
```

If your batch directory follows the default layout, this shorter command is equivalent:

```bash
python scripts/run_remove_slipknots.py -b batch_results
```

That command reads clustered CSVs from `batch_results/clustered`, writes cleaned CSVs to `batch_results/no_slipknot`, and writes a summary file:

```text
batch_results/no_slipknot/remove_slipknots_summary.csv
```

## End-to-end example

This example uses `2ww4.pdb` from the repository's `test/PDB` directory. Clone the repository first to obtain that test structure:

```bash
git clone https://github.com/vuqv/glink_entanglement.git
cd glink_entanglement
pip install -e .

mkdir -p out clustered

glink \
  -f test/PDB/2ww4.pdb \
  -o out \
  --GLN_threshold 0.6 \
  --contact_cutoff 4.5 \
  --topoly_density 0 \
  --topoly_min_dist 10 6 5

glink-cluster \
  -f out/2ww4_glink_contacts.csv \
  -o clustered

python scripts/remove_slipknots.py \
  -i clustered/2ww4_glink_contacts_clustered.csv \
  -o cleaned
```

Inspect the outputs:

```bash
head out/2ww4_glink_contacts.csv
head clustered/2ww4_glink_contacts_clustered.csv
head cleaned/2ww4_glink_contacts_clustered_no_slipknot.csv
```

Count rows:

```bash
wc -l out/2ww4_glink_contacts.csv
wc -l clustered/2ww4_glink_contacts_clustered.csv
wc -l cleaned/2ww4_glink_contacts_clustered_no_slipknot.csv
```

Remember that `wc -l` includes the header.

## Batch processing

For a directory of PDB files, the batch wrapper can run the raw `glink` step and the clustering step together:

```bash
python scripts/run_batch.py \
  -i PATH_TO_PDB_DIR \
  -o batch_results \
  -j 8
```

The batch output layout is:

```text
batch_results/
├── raw/
├── clustered/
└── batch_summary.csv
```

For each input PDB, the script first writes:

```text
batch_results/raw/PDB_STEM_glink_contacts.csv
```

Then it checks whether the raw CSV has at least one result row. If raw entanglements are present, it runs `glink-cluster` and writes:

```text
batch_results/clustered/PDB_STEM_glink_contacts_clustered.csv
```

If the raw CSV has only a header and no entanglement rows, the batch summary records `status` as `no_entanglement` and does not run clustering for that PDB.

The wrapper also accepts a single PDB file:

```bash
python scripts/run_batch.py \
  -i PATH_TO_PDB \
  -o batch_results
```

Use `--force` to recompute outputs that already exist:

```bash
python scripts/run_batch.py \
  -i PATH_TO_PDB_DIR \
  -o batch_results \
  --force
```

The default number of workers is the available CPU count. Set `-j` when you want to limit concurrency for a workstation or job allocation.

After batch clustering, optional slipknot cleanup can be run over the whole batch:

```bash
python scripts/run_remove_slipknots.py -b batch_results
```

If your inputs are keyed by UniProt ID, the repository also includes a table-driven wrapper. It reads a CSV, TSV, pickle DataFrame, or similar table with a `Uniprot` column by default; each value is resolved to `<pdb_dir>/<Uniprot>.pdb`.

```bash
python scripts/run_batch_uniprot.py \
  -t PROTEIN_TABLE.csv \
  -p PATH_TO_PDB_DIR \
  -o batch_results \
  -j 8
```

Use `--column` if the UniProt column has a different name, `--delimiter` for nonstandard text delimiters, and `--extension` if the structure files do not end in `.pdb`.

You can still run a manual shell loop when you want direct control over every command:

```bash
mkdir -p out clustered

for pdb in PATH_TO_PDB_DIR/*.pdb; do
  glink -f "$pdb" -o out
done

for csv in out/*_glink_contacts.csv; do
  if [ "$(wc -l < "$csv")" -gt 1 ]; then
    glink-cluster -f "$csv" -o clustered
  fi
done
```

## Using the Python API

The command-line tools format their output for convenient CSV use. The underlying functions can also be imported into larger analysis pipelines.

To calculate raw contacts:

```python
from glink_entanglement.glink import calculate_pdb_glink

df = calculate_pdb_glink(
    "PATH_TO_PDB",
    threshold=0.6,
    cutoff=4.5,
    topoly_density=0,
    topoly_min_dist=(10, 6, 5),
)

print(df.head())
```

The Python-level DataFrame returned by `calculate_pdb_glink` includes internal residue fields such as `chain`, `resid_i`, `resname_i`, `resid_j`, `resname_j`, `contact_i_index`, and `contact_j_index`. The CLI converts those fields into the public `contact`, `i`, and `j` columns before writing CSV.

To create the same public-style CSV from Python:

```python
output = df.copy()
output.insert(
    0,
    "contact",
    output["resname_i"] + output["resid_i"].astype(str) + "-" + output["resname_j"] + output["resid_j"].astype(str),
)
output = output.drop(columns=["chain", "resid_i", "resname_i", "resid_j", "resname_j"])
output = output.rename(columns={"contact_i_index": "i", "contact_j_index": "j"})
output.to_csv("out/PDB_STEM_glink_contacts.csv", index=False)
```

To cluster a raw CSV:

```python
from glink_entanglement.clustering import cluster_glink

clustered = cluster_glink(
    "out/PDB_STEM_glink_contacts.csv",
    cutoff=52,
)

clustered.to_csv("clustered/PDB_STEM_glink_contacts_clustered.csv", index=False)
```

## Parameter guidance

Start with the defaults, then change one parameter at a time.

Use `--contact_cutoff` to change the native-contact definition. Larger cutoffs produce more candidate contacts and may increase runtime. Smaller cutoffs are stricter and may miss looser loop-closing contacts.

Use `--GLN_threshold` to change how easily partial Gaussian-linking values become nonzero `Gn` and `Gc` scores. A lower threshold is more permissive; a higher threshold is stricter.

Use `--topoly_density 1` when you need Topoly's denser default surface calculation and can tolerate the slower run.

Use `--topoly_min_dist` when you need to tune Topoly's crossing-reduction behavior. The three values are crossing, loop, and tail-end distances.

Use `glink-cluster --cutoff` when the built-in organism presets are not appropriate for the structures being compared.

## Interpreting results

A row in the raw CSV means:

- two residues form a heavy-atom contact;
- the contact closes a loop with nonzero rounded Gaussian linking against at least one side of the chain;
- Topoly found at least one crossing residue on a side with nonzero rounded linking.

A row in the clustered CSV means:

- one or more raw contacts were judged to represent the same entanglement signature;
- the reported `contact` is the representative loop-closing contact;
- `num_contacts` gives the number of raw contacts summarized by that representative;
- `contacts` lists the raw contact labels that contributed to the representative.

For downstream comparisons across proteins or models, start with the clustered CSV. For method debugging, parameter tuning, or manual inspection of a specific entanglement, return to the raw CSV.

## Troubleshooting

If `glink` writes zero contacts, check the structure before changing many parameters:

- Are protein atoms recognized by MDAnalysis?
- Are C-alpha atoms present for the analyzed residues?
- Are chain breaks, insertion codes, or missing residues expected?
- Is the contact cutoff too small?
- Is the GLN threshold too strict?
- Did Topoly fail to find crossings on the sides with nonzero `Gn` or `Gc`?

If all atoms appear under chain `SYSTEM`, the PDB may not contain explicit chain or segment identifiers. The calculation can still run, but adding chain identifiers helps interpretation when a file contains multiple chains.

If `glink-cluster` reports missing columns, make sure its input is the raw CSV produced by `glink`, not a clustered CSV or a manually edited table. The current raw public schema should include:

```text
contact,i,j,gn,gc,Gn,Gc,crossingsN,crossingsC,crossings
```

If the shell cannot find `glink` or `glink-cluster`, confirm that the active shell is using the Python environment where the package was installed:

```bash
python -m pip show glink-entanglement
which glink
which glink-cluster
```

After reinstalling in the same environment, refresh the shell's command cache:

```bash
hash -r
```

## Summary

`glink-entanglement` turns one all-atom protein PDB into two practical tables: a raw contact-level entanglement table from `glink`, and a representative entanglement table from `glink-cluster`. Optional scripts can then remove slipknot crossing pairs or batch the workflow across many PDB files. The current output contract is centered on readable contact labels such as `ALA10-GLY42`, zero-based C-alpha contact indices `i` and `j`, side-specific Topoly crossing columns `crossingsN` and `crossingsC`, and semicolon-delimited raw-contact membership in the clustered output.
