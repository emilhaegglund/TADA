# SDBW
This is a Snakemake workflow to subsample biological sequence database based on taxonomy and create blast database from them.

## Installing
To run the workflow, first clone it from git using
```
git clone https://github.com/emilhaegglund/SDBW.git
```

## Setting up the configuration file
The first thing to do before running the workflow is to setup the
configuration-file. This is located in the `config/config.yaml`.
This determines how the workflow will behave.
The first option is to set the path to the output-directory:
```
path:
  results: "results"
```
There are three main different methods; first is to use the NCBI Taxonomy to subsample the RefSeq database, the second method is to sample based on the GTDB Taxonomy, the final method is to do the sampling to retain the largest evolutionary distances between the genomes, this is done by trimming the phylogenies from GTDB. The main differences is that RefSeq will also include Eukaryotes, while GTDB is limited to bacteria and archaea. Sampling GTDB can also be combined with filtering on a
large number of metadata associated with each entry.
To determine the methods to use
```
method:
  sample_refseq: False
  sample_gtdb: True
  trim_gtdb: True
```
### Options for sampling RefSeq

### Options for sampling GTDB
Similar to the RefSeq, the main option is to set the taxonomical level to sample at, and also the maximum number of taxa to include.
```
sample_gtdb:
  max_taxa: 3
  taxonomic_level: "family"
```
The above example will sample 3 entries from each family.
The GTDB-database comes with a large collection of metadata associated with each taxa. This can be use to first exclude entries. Currently

TODO:
### Building a sampling scheme
We can also do a more advanced sampling by adding the path to a file containing a sampling scheme.

```
phylum: Name
  max_taxa: Number
```

### Option for trimming GTDB
There is two differet phylogenies, one for bacteria and one for archaea. At the moment trim the phylogenies until each of them contains the number of taxa as defined below.
__Alternative method:__ Only define the number of bacteria, and start trimming this phylogeny. Save the evolutionary distance for the last trimmed leaf-pair and then trim the Archaeal tree until we reach the same phylogenetic distance.
```
trim_gtdb:
  number_bac: 1000
  number_ar: 200
```


## Running the workflow
