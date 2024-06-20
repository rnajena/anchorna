### AnchoRNA
[![build status](https://github.com/rnajena/anchorna/workflows/tests/badge.svg)](https://github.com/rnajena/anchorna/actions)
[![codecov](https://codecov.io/gh/rnajena/anchorna/branch/master/graph/badge.svg)](https://codecov.io/gh/rnajena/anchorna)

Find anchors in RNA/DNA sequences.

Anchorna is designed to find anchor regions in a set of similar sequences. Anchors are saved as GFF files.
Anchorna provides export functionality to Jalview feature files. You can cut out sequences relative to anchors for import into external tools.

Anchorna supports curation and iterative refinement of anchors.

### Installation

Use pip to install anchorna. It is best installed in a dedicated (conda) environment.
To install the latest release use the following command and replace `VERSION` with the latest version number.

```
pip install https://github.com/rnajena/anchorna/archive/refs/tags/VERSION.zip

```

To install the development version in editable mode, clone the repository and install with `pip -e repo_path`.

### Usage & Tutorial

A command line tool is provided, check `anchorna -h` and the help of subcommands. An example configuration file can be created with `anchorna create`.

Start the tutorial with

```
anchorna create --tutorial
```

We provide here some CLI calls. Please run each of them separately and carefully try to understand what they do.

```
anchorna go anchors.gff

# Print anchors
anchorna print anchors.gff
anchorna print -v anchors.gff
anchorna print --mode seq anchors.gff

# View anchors in Jalview (amino acid sequences)
anchorna export anchors.gff -o jalview_fts.txt
sugar translate --cds pesti_example.gff -o pestiaa.fasta
jalview pestiaa.fasta --features jalview_fts.txt
# Instead of the above three lines we can just use anchorna view
anchorna view anchors.gff

# View anchors in Jalview (nucleotide sequences)
anchorna view --mode seq anchors.gff

# We are not satisfied with the first anchor, remove it and check the resulting anchors in Jalview, afterwards save to new anchor file
anchorna combine "anchors.gff||a0" | anchorna view -
anchorna combine "anchors.gff||a0" -o anchors2.gff
# Note that anchors are renumbered, the anchor number is not a property of each anchor,
# it is just used as a simple reference later on
anchorna print anchors2.gff

# Cutout subsequences starting before start codon and
# ending 10 nucleotides after 3rd anchor (remember 0-based indexing) for usage in external tool
anchorna cutout anchors2.gff "ATG<" "A2>+10" -o pestiA3.fasta


### Advanced
# Load anchors into IPython session
anchorna load anchors2.gff

# Run anchorna on subsequences and combine results into a single anchor file
# We are not satisfied with the long region (>1000 aa) without any anchors between anchors A52 and A53
anchorna cutout anchors2.gff "A52>" "A53<" -o pesti_sub.sjson  # Need to use sjson format for now, to save the offsets
anchorna go --fname pesti_sub.sjson --thr-quota 0.9 anchors_sub.gff  # adapt options
# We found another 38 anchors and investigate these
anchorna print anchors_sub.gff
anchorna view anchors_sub.gff
# Finally we merge the two anchor files
anchorna combine anchors2.gff anchors_sub.gff -o anchors_final.gff
anchorna print anchors_final.gff
anchorna view anchors_final.gff
```

### Run tests

Clone or download this repository cd into anchorna folder and call `pytest`.
