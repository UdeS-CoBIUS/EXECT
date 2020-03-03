# EXECT
Program Version 2020

Reconstructs the Phylogeny of a set of bacterial genomes while accounting for horizontal gene transfers.
----------------------------------------------------------------

Authors: Nilson Da Rocha Coimbra, and Aida Ouangraoua

Université de Sherbrooke, Canada
CoBIUS Lab:  https://cobius.usherbrooke.ca/

For questions email us at aida.ouangraoua@usherbrooke.ca

### Requirements:

- EXECT has been tested on Linux and OSX

- Python 3.7 or higher  

- Python modules: os, glob, time, decimal, ete3, pandas, matplotlib, seaborn, Bio

- Mafft+ standalone

- IslandPath-Dimob

- Orthofinder

- FastTree

- RAxML

- ASTRAL

- ASTRID

- CONSENSE (included in PHYLIP package)

### Installation

The easiest way to run EXECT on Linux or OSX is to use Bioconda, and
install IslandPath, OrthoFinder, FastTree, RAxML and PHYLIP via bioconda

exemple : conda install islandpath orthofinder fasttree raxml phylip

For IslandPath, a release is included in the repository in the
`tools/islandpath` directory. You can download the most recent
version at https://github.com/brinkmanlab/islandpath/releases/

For ASTRAL, version 5.5.9 is included in the repository in the
`tools/ASTRAL-master` directory. You can download the most recent
version at https://github.com/smirarab/ASTRAL/

For ASTRID : Download the appropriate binary file for your operating
system (ASTRID-osx or ASTRID-linux) at
https://github.com/pranjalv123/ASTRID/releases, and put it in your path.


### Usage

#### First, edit the file `Utils.py` to update the values of variables:

DATADIR : directory where the results will be written

EMAIL : email address for connexion to NCBI databases

GENOMELISTFILE : list of genomes to include in the analysis

ISLANDPATHDIR : path to IslandPath release

ASTRALJARFILE : path to ASTRAL jar file

ASTRIDEXECUTABLENAME : Name of ASTRID binary file

#### Then, run the program
```
python3.7 ComputeSpeciesTree.py

```

