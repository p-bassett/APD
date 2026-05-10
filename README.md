# Analysing amyloid fibril structures for their sidechain packing interactions

>This is a modified fork of [APD](https://github.com/3dem/APD.git) published by S.H.W. Scheres. This fork is designed with modifications for batch-processing of structures and creating figures analyzing the distance in structural similarity measured by APD.

The Amyloid Packing Difference (APD) is a metric for comparison of side-chain packing interactions in amyloid structures. Calculating the APD is done by first mapping contacts between residues in a single structure, then comparing the contacts between separate structures in a pairwise manner. More information about the analysis can be found in the corresponding preprint from S.H.W. Scheres [here](https://doi.org/10.64898/2026.02.18.706523).

## Purpose of this project
This repository contains additional code as well as modified versions of the original scripts to allow analyzing several amyloid fibril structures in a high-throughput manner. The analysis then allows for performing hierarchical clustering based on APD value as the distance metric.

## Modifications to the workflow
- Both `compare.py` and `contacts.py` were restructured in order to be called more easily in a new script, `batch-apd.py`, which performs the APD analysis end-to-end on each PDB structure in an input directory `PDBs/`.
- The output of `batch-apd.py` varies slightly in structure from the original scripts, automatically grouping similar files together (such as .log and .svg files).
- The `batch-apd.py` script creates two distance matrices using the XY and XYZ APD values for each pairwise structural comparison. These can be used as-is or for downstream analysis.
- An additional `cluster.py` script takes the distance matrices from `batch-apd.py` and performs hierarchical clustering to sort structures by their APD values. A clustering dendrogram and a heatmap are both generated.
- The underlying logic of the APD calculations is unchanged.

## Installation
Necessary dependencies are **Biopython, Matplotlib, SciPy,** and **Seaborn**. These can be installed manually or via Conda using the environments.yml file.

```bash
# Install using pip
pip install biopython matplotlib scipy seaborn

# Install using conda
conda env create -f environments.yml
```

## Running the scripts
The original `compare.py` and `contacts.py` scripts can be used via the command line using the same original usage.

```bash
# Determine contacts of a structure
python contacts.py <InputPDB.pdb>

# Compare contacts between two structures
python compare.py <PDB1contacts.csv> <PDB2contacts.csv>
```

The batch-processing script automatically detects any structures found in `PDBs/` and iterates through each of them.

```bash
python batch-apd.py
```

The clustering script automatically finds the distance matrices from `batch-apd.py` in the output `Comparisons/` directory.

```bash
python cluster.py
```

## Fine-tuning figures
The original `compare.py` script, when used via the command line, allows for reorienting the compared structures in the .svg figure using flags such as `--rot1` or `--flip2`. See the original repo for more details.

When using the `batch-apd.py` script, all structures are aligned to a common orientation to keep the pairwise comparison figures as consistent as possible. This orientation for all structures can be adjusted using similar flags:

```bash
python batch-apd.py --rotate_all <degrees> --flip_all
```

The `compare.py` script has been modified to generate residue numbering to mark specific residue positions.

The `compare.py` script has been modified to generate residue numbering to mark specific residue positions.

## References
Scheres, S. H. W. The amyloid packing difference: a pairwise comparison metric for amyloid structures. _bioRxiv_ 2026.02.18.706523 Preprint at https://doi.org/10.64898/2026.02.18.706523 (2026).

## License
This modified fork is distributed under the MIT License, in accordance with the original APD repository authored by Sjors H.W. Scheres (MRC Laboratory of Molecular Biology). See the `LICENSE` file for more details.