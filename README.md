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

### Defining a sampling scheme
The subsampling of Refseq or GTDB is based on a sampling scheme that is defined in a YAML-file with the following structure. Some examples of sampling schemes can be found in the config-directory.
```
taxonomic_name:
    sampling_level: [domain, phylum, class, order, family, genus, species]
    taxa: <int>
```
The `taxonomic_name` represents the name of the taxa we would like to sample from. It has to be a name that is defined in the GTDB or NCBI taxonomy, depending on which database to sample from. The key-word `all` can be used to sample from all Bacteria and Archaea. The `sampling_level` should be followed by a valid taxonomic level, this is the level at which the sampling will be done at. Finally, the `taxa` should be followed by an integer defining the number of taxa we want to sample from each group.
For example, the sampling scheme below will sample 10 taxa from each phylum in both Bacteria and Archaea.
```
all:
    sampling_level: phylum
    taxa: 10
```
We can also create more advanced sampling scheme by inluding several sampling criteria. This will do a hierarchical sampling starting sampling at the lowest taxonomic level and continue to sample at higher and higher levels.
Taxonomic groups that have already been sampled from is excluded from sampling at higher levels.

In the example below we will sample 10 taxa from `Bartonella bacilliformis`, for remaining species in the `Bartonella` genus we will sample 2 taxa per species. Finally, to this we will sample 5 taxa from the `Rhizobiales_A` order outside of the `Bartonella` genus.
```
Bartonella bacilliformis:
  sampling_level: species
  taxa: 10

Bartonella:
  sampling_level: species
  taxa: 2

Rhizobiales_A:
  sampling_level: order
  taxa: 5
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
GTDB uses the latest release r207. Before any subsampling is done, records that have been suppressed from RefSeq falls back to the corresponding GenBank record. Records that have been suppressed from GenBank is removed from the metadata-files from GTDB and will be excluded from the sampling.
The GTDB-database comes with a large collection of metadata associated with each taxa. This can be use to first exclude entries. Currently, the subsampling- and pruning-workflow support exclusion of taxa that doesn't fullfill completeness and contamination criteria defined in the config file. The sample-workflow also support the exclusion of non-gtdb-representative taxa.

### Options for subsampling of GTDB
```
subsample_gtdb:
    sampling_scheme: <path>
    completeness: <int>
    contamination: <int>
    gtdb_species_representative: <bool>

```
Options:
    `sampling_scheme:` path to the sampling scheme that will be used.
    `completeness:` An integer between 0-100, will exclude entries with a completeness estimate less than this value.
    `contamination:` An integer, will exclude entries with a contamination estimate larger than this value.
    `gtdb_species_representative:` False will keep all entries while True will only keep entries that are classified as GTDB species representatives.

In the example below we will use the sampling scheme defined in `config/sampling_scheme.basic.yaml`. First, taxa with an estimated completeness below 90%, an estimated contamination over 5%, and taxa that is not gtdb-representative will be excluded. From the remaining taxa we will sample one representative for each phylum defined in GTDB.
```
subsample_gtdb:
  sampling_scheme "../config/sampling_scheme.basic.yaml"
  completeness: 90
  contamination: 5
  gtdb_species_representative: True
```

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
