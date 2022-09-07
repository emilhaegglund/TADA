# SDBW
A snakemake workflow to subsample taxa from GTDB or RefSeq based on taxonomic/phylogenomic information and create a protein-sequence blast database (Diamond) from them.

## Installing
First clone this repository from git using
```
git clone https://github.com/emilhaegglund/SDBW.git
```
Next install and activate a conda environment, this will install mamba and snakemake.
```
conda env create -f environment.yml
conda activate sdbw
```

## Setting up the configuration file
The first thing to do before running the workflow is to setup the
configuration-file that will determine how the workflow will behave. An example is located in the `config/config.yaml`.
You can either modify this file or create a new, but then you will have to
specify the location in the `--configfile` option when running snakemake.

The first option is to set the path to the output-directory:
```
path:
  results: "results"
```
There are three main different methods:
1. Use the NCBI Taxonomy to subsample the RefSeq database.
2. Subsample based on the GTDB Taxonomy.
3. Subsampling to retain the largest evolutionary distances between the taxa. This is done by pruning the phylogenies from GTDB.

The main differences is that RefSeq will also include Eukaryotes, while GTDB is limited to bacteria and archaea. Sampling GTDB can also be combined with filtering on a large number of metadata associated with each entry.
To determine the methods to use change the following lines in the config-file.
```
method:
  subsample_refseq: False
  subsample_gtdb: True
  prune_gtdb: False
```
### Options for sampling RefSeq
__TODO__: Not implemented yet.

### GTDB
GTDB uses the latest release r207. Before any subsampling is done, we remove records that have been suppressed from GenBank.

### Options for subsampling of GTDB
Similar to the RefSeq, the main option is to set the taxonomical level to perform the subsampling, and also the maximum number of taxa to include. Valid taxonomical levels are `domain,phylum,class,order,family,genus,species`.
```
subsample_gtdb:
  max_taxa: 3
  taxonomic_level: "family"
```
The above example will sample 3 entries from each family defined in GTDB.

The GTDB-database comes with a large collection of metadata associated with each taxa. This can be use to first exclude entries. Currently, the subsampling- and pruning-workflow support exclusion of taxa that doesn't fullfill completeness and contamination criteria defined in the config file. The sample-workflow also support the exclusion of non-gtdb-representative taxa.
```
subsample_gtdb:
  max_taxa: 1
  taxonomic_level: "phylum"
  completeness: 90
  contamination: 5
  gtdb_representative: True
```
The example above will first remove taxa with an estimated completeness below 90%, an estimated contamination over 5%, and taxa that is not gtdb-representative.

### Option for pruning the GTDB-phylogeny
GTDB contains one phylogeny for bacteria and one for archaea. SDBW will prune the phylogenies based on the phylogenetic distance between taxa, until the phylogeny contains the number of taxa as defined in the config file.
Before the distance-based pruning we can also use the `completeness` and `contamination` options to first remove taxa that does not fullfill this criteria.

__Alternative method 1:__ Only define the number of bacteria, and start trimming this phylogeny. Save the evolutionary distance for the last trimmed leaf-pair and then trim the Archaeal tree until we reach the same phylogenetic distance.

__Alternative method 2:__ Just give the total number of taxa to retain, then trim both the bacteria and archaea phylogny in parallel, can be a bit more technical to achieve.
```
prune_gtdb:
  bac120: 1000
  ar53: 200
  completeness: 90
  contamination: 5
```
The above example will prune all taxa that has a completeness of less than 90% and a contamination over 5%, then it will continue to prune the bacterial phylogeny untill 1000 taxa remains, and the archaeal phylogeny untill 200 taxa remains.

## Running the workflow
Finally, to run the workflow, use the following command. If you have the config-file in a
different location, replace the path after `--configfile`
```
cd workflow
snakemake --cores 4 --use-conda --conda-frontend mamba --configfile ../config/config.yaml
```
