# TADA - Taxonomic Aware Dataset Assembly
A Snakemake workflow to assemble dataset for comperative and phylogenetic analysis of bacteria and archaea. Datasets can be generated baed on taxonomic information from [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy), or taxonomic or phylogenomic information from [GTDB](https://gtdb.ecogenomic.org).

## Installing
Running the TADA-workflow requires Conda. First clone the repository from git
and change into the TADA directory.

```
git clone https://github.com/emilhaegglund/TADA.git
cd TADA
```

Next, install and activate the conda environment from which the
workflow will be run. This will install Mamba and Snakemake.

```
conda env create -f environment.yaml
conda activate tada
```

## Setting up the configuration file
Before running the workflow, the first step is to set up the configuration file. This file will determine the behavior of the workflow. An example of the configuration file can be found in `config/config.yaml`. You can either modify this file or create a new. The location of the config-file must be specified using the `--configfile` option when running Snakemake.

The first option is to set the path to the output-directory:
```
workdir: "results"
```

### Sampling methods
The workflow can be run using three different methods:
* Sampling based on the NCBI taxonomy.
* Sampling based on the GTDB taxonomy.
* Sampling based on the GTDB phylogeny.

To determine which sampling method to use select one of the following option: `sample_ncbi`, `sample_gtdb`, or `prune_gtdb`.

```
method: "sample_gtdb"
```

A random seed can be used to reproduce the output of sampling and pruning from the GTDB-database.

```
seed:42
```

### Select output
There are several output options for TADA. First we can select if the workflow should download genomes, cds, or proteomes for the sampled genomes. If all options below are set to False, the workflow will stop after the sampling procedure. If an annotation does not exist for a sampled taxa, the genome will be annotated using Prokka.

```
download_genomes: False
download_cds: False
download_proteomes: True
```

We can also make the workflow build different type of Blast-databases, either using the [NCBI Blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables) suite or [Diamond](https://github.com/bbuchfink/diamond).

```
output:
    blast_genomes: False
    blast_cds: False
    blast_protein: False
    diamond_protein: True
```

### Options for sampling
Next follows option that are specific to the different sampling methods listed above.

__Options for sampling the GTDB Taxonomy__
```
sample_gtdb:
    sampling_scheme: <path>
    completeness: <float>
    contamination: <float>
    gtdb_species_representative: <bool>
    version: <str>
```

`sampling_scheme`: Path to the sampling scheme that will be used. See Defining a sampling scheme for more details on this.

`completeness`: Exclude taxa with a completeness estimate less than this value (Default: `0`).

`contamination`: Exclude taxa with a contamination estimate larger than this value (Default: `100`).

`gtdb_species_representative`: False will keep all entries while True will only keep entries that are classified as GTDB species representatives.

`version:` Select which version of GTDB to use, `207` and `214` are supported (Default: `214`).

__Example__

In the example below we will use the sampling scheme defined in `config/sampling_scheme.basic.yaml`. For this example the workflow will sample three taxa for each phylum, but only sample from representative species with an estimated completeness over 90% and an estimated contamination under 5%.

```
subsample_gtdb:
  sampling_scheme "../config/sampling_scheme.basic.yaml"
  completeness: 90
  contamination: 5
  gtdb_species_representative: True
```

__Options for pruning the GTDB phylogenies__

GTDB includes separate phylogenies for bacteria and archaea. TADA will prune these phylogenies based on the evolutionary distance between taxa, reducing the number of taxa to the amount specified in the configuration file. Before the distance-based pruning, it is also possible to use the completeness and contamination criteria to remove taxa that do not meet these requirements from the phylogeny.

```
prune_gtdb:
    bac120: <int>
    ar53: <int>
    completeness: <float>
    contamination: <float>
    prune_method: <str>
    version: <str>
```

`bac120`: Number of taxa to sample from the bacterial phylogeny.

`ar53`: Number of taxa to sample from the archaeal phylogeny.

`completeness`: Exclude taxa with a completeness estimate less than this value (Default: `0`).

`contamination`: Exclude taxa with a contamination estimate larger than this value (Default: `100`).

`prune_method`: Select what method to use for pruning, `"shortest"` will keep the taxa with the shortest branch in a leaf-pair, `"longest"` will keep the taxa with the longest branch in a leaf-pair, and `"random"` will randomly select the taxa to keep in a leaf-pair (Default: `"shortest"`).

`version:` Select which version of GTDB to use, `207` and `214` are supported (Default: `214`).

__Example__

In the example belwo TADA will first remove all taxa with an estimated completeness under 90% and an estimated contamination over 5%. It will then continue to prune the bacterial phylogeny untill 1000 taxa remains. For the archaeal phylogeny it will prune the phylogeny until 200 taxa remains.

```
prune_gtdb:
  bac120: 1000
  ar53: 200
  completeness: 90
  contamination: 5
  prune_method: "shortest"
```

__Options for sampling from the NCBI Taxonomy__

To sample from the NCBI Taxonomy we have to give the path to a sampling scheme and we also need to define if we want to sample from GenBank or RefSeq. Sampling from NCBI is restricted to taxa classified as Bacteria or Archaea. The reason for this is that the annotation software in the workflow are for prokaryotic genomes.
```
sample_ncbi:
    sampling_scheme: <path>
    source: <string>
```
`sampling_scheme`: Path to the sampling scheme that will be used. See Defining a sampling scheme for more details on this.

`source`: Sample from `"GenBank"` or `"RefSeq"`.

__Example__

In the example below we will sample one taxa from each defined phylum in the RefSeq-database.

```
sample_ncbi:
    sampling_scheme: "../config/sampling_scheme.ncbi_refseq.yaml"
    source: RefSeq
```

## Defining a sampling scheme

The sampling is based on a sampling scheme that is defined in a YAML-file with a structure described below. Examples of sampling schemes can be found in the config-directory.

```
taxonomic_name:
    sampling_level: [domain, phylum, class, order, family, genus, species]
    taxa: <int> or "all"
```
`taxonomic_name`: The name of the taxa to sample from. The name has to be  defined in the GTDB or NCBI taxonomy, depending on which database to sample from. The key-word `all` can be used to sample from all Bacteria and Archaea.

`sampling_level`: The taxonomic level to perform the sampling at.

`taxa`: The number of taxa we want to sample from each group at that taxonomic level. The key-word `all` can also be used, this will keep all taxa in the groups.

__Example: Basic sampling scheme__

The sampling scheme below will sample 10 taxa from each phylum in both Bacteria and Archaea.
```
all:
    sampling_level: phylum
    taxa: 10
```

__Example: Complex sampling scheme__

It is also possible to construct more complex sampling scheme by including multiple sampling criterias. The workflow will then perform a hierarchical sampling starting the sampling procedure at the lowest taxonomic level and then continue  to sample at higher and higher levels. Taxonomic groups which have already been used to sample from is excluded from sampling at higher levels.

The example below demonstrates a complex sampling scheme. We will begin by sampling 10 taxa from the species Bartonella bacilliformis. For the remaining species within the Bartonella genus we will sample two taxa per species. Adding to this, we will sample five taxa from the Rhizobiales_A order outside of the Bartonella genus.

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

## Running the workflow
To run the workflow, use the following command. If you have the config-file in a different location, replace the path after --configfile

```
cd workflow
snakemake --cores 4 --use-conda --conda-frontend mamba --configfile ../config/config.yaml
```
