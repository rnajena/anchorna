### AnchoRNA
[![build status](https://github.com/rnajena/anchorna/workflows/tests/badge.svg)](https://github.com/rnajena/anchorna/actions)
[![codecov](https://codecov.io/gh/rnajena/anchorna/branch/master/graph/badge.svg)](https://codecov.io/gh/rnajena/anchorna)
[![docs status](https://readthedocs.org/projects/anchorna/badge/?version=latest)](https://anchorna.readthedocs.io)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12190267.svg)](https://doi.org/10.5281/zenodo.12190267)

Find anchors in RNA/DNA sequences.

AnchoRNA is designed to find anchors, i.e. conserved regions, in coding sequences in a set of similar sequences. Anchors are stored as GFF files.
Anchors can be viewed in JalView. You can cut out sequences relative to anchors for import into external tools.

AnchoRNA supports curation and iterative refinement of anchors. It can also be used to find anchors in non-coding sequences or amino acid sequences.

### Installation

Use pip to install anchorna. It is best installed in a dedicated (conda) environment.
To install the latest release use the following command.

```
pip install https://github.com/rnajena/anchorna/archive/refs/tags/1.0.1.zip
```

To install the latest development version, use:

```
pip install https://github.com/rnajena/anchorna/archive/refs/heads/master.zip
```

To install the package in editable mode, clone the repository and install with `pip -e repo_path`.

### Usage & Tutorial

A command line tool is provided, see `anchorna -h` and the help of subcommands. An example configuration file can be created with `anchorna create`.
See the example configuration file for a description of AnchoRNA's options. Consult the (API documentation)[https://anchorna.readthedocs.io] for details.

Start the tutorial with

```
anchorna create --tutorial
```

We provide here some CLI calls. Please run each of them separately and carefully try to understand what they do.

```
# Identify anchors, use 8 cores
anchorna go --njobs 8 anchors.gff

# Print anchors
anchorna print anchors.gff
anchorna print -v anchors.gff
anchorna print --mode nt anchors.gff

# View anchors in Jalview (amino acid sequences) and
# align with the second anchor
anchorna view anchors.gff --align A1
# For reference, the following commands are called under the hood
# by anchorna view
anchorna export anchors.gff --fmt jalview -o jalview_fts.txt
sugar translate --cds pesti_example.gff -o pestiaa.fasta
jalview pestiaa.fasta --features jalview_fts.txt

# View anchors in Jalview (nucleotide sequences)
anchorna view --mode nt anchors.gff

# We are not satisfied with the first anchor, remove it and check the resulting anchors in
# Jalview, afterwards save to new anchor file
anchorna combine "anchors.gff||a0" | anchorna view -
anchorna combine "anchors.gff||a0" -o anchors2.gff
# Note that anchors are renumbered, the anchor number is not a property of each anchor,
# it is just used as a simple reference later on
anchorna print anchors2.gff

# Cutout subsequences starting before start codon and ending 10 nucleotides
# after 3rd anchor (remember 0-based indexing) for usage in external tool
anchorna cutout anchors2.gff "ATG<" "A2>+10" -o pestiA3.fasta


### Advanced
# Load anchors into IPython session
anchorna load anchors2.gff

# Run anchorna on subsequences and combine results into a single anchor file
# We are not satisfied with the long region (>1000 aa) without any anchors,
# find two anchors immediately before and after this "empty" region,
# replace the "??" signs with the anchor numbers
anchorna view anchors2.gff
anchorna cutout anchors2.gff "A??>" "A??<" -o pesti_sub.sjson  # Need to use sjson format for now, to save the offsets
anchorna go --fname pesti_sub.sjson --thr-quota-add-anchor 0.9 anchors_sub.gff  # adapt options
# We found more anchors and investigate these
anchorna print anchors_sub.gff
anchorna view anchors_sub.gff
# Finally we merge the two anchor files
anchorna combine anchors2.gff anchors_sub.gff -o anchors_final.gff
anchorna print anchors_final.gff
anchorna view anchors_final.gff
```

### Run tests

Clone or download this repository cd into tests folder and call `pytest`.
