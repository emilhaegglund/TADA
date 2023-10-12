# TADA - Taxonomic Aware Dataset Assembly
A Snakemake workflow to assemble balanced, representative and manageable datasets for comparative and phylogenetic analysis of bacteria and archaea. Datasets can be generated based either on the phylogenomic tree offered by [GTDB](https://gtdb.ecogenomic.org) , or on the taxonomy offered by GTDB or by [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy).

## Dependency
Running the TADA-workflow requires [Conda](https://docs.conda.io/projects/conda/en/stable/).


## Installing
Clone the repository from git and change into the TADA directory.

```
git clone https://github.com/emilhaegglund/TADA.git
cd TADA
```

Install and activate the conda environment from which the
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

### Choice of sampling method
The workflow can be run using three different methods:
* Sampling based on the NCBI taxonomy (`sample_ncbi`).
* Sampling based on the GTDB taxonomy (`sample_gtdb`).
* Sampling based on the GTDB phylogeny (`prune_gtdb`).

E.g.:

```
method: "sample_gtdb"
```

A random seed can be used to reproduce the output of sampling and pruning from the GTDB-database.

```
seed: 42
```

When using the `sample_gtdb` or `sample_ncbi` option a file containing a list of genome accessions to be include in the dataset can be given with the `required` option.
```
required: "../config/required-genomes.txt"
```

An NCBI API Key can be used in the workflow with the following option
```
ncbi_api_key: "NCBI-API-KEY"
```

### Select what to download
TADA can download genomes, CDS (genes), proteomes, and/or GFF3 annotations for the sampled genomes. If all options below are set to False, the workflow will stop after the sampling procedure. TADA will annotate genomes for which no annotation is available using [Prokka](https://github.com/tseemann/prokka).

```
downloads:
  genomes: False
  cds: False
  proteomes: True
  gff3: False
```

### Select what databases to create
TADA can also build different type of Blast-compatible databases, either using the [NCBI Blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables) suite or [Diamond](https://github.com/bbuchfink/diamond) (only for proteins).

```
databases:
    blast_genome: False
    blast_cds: False
    blast_protein: False
    diamond_protein: True
```

### Options for sampling
Next follows options specific to the different sampling methods listed above.

__Options for sampling from the NCBI Taxonomy__

To sample from the NCBI Taxonomy we have to give the path to a sampling scheme and we also need to define if we want to sample from GenBank or RefSeq. Sampling from NCBI is restricted to taxa classified as Bacteria or Archaea. The reason for this is that the annotation software in the workflow are for prokaryotic genomes.
```
sample_ncbi:
    sampling_scheme: <path>
    database: <string>
```
`sampling_scheme`: Path to the sampling scheme that will be used. See [Defining a sampling scheme](#defining-a-sampling-scheme) for more details on this.

`database`: Sample from `"GenBank"` or `"RefSeq"`.

__Example__

In the example below, TADA will sample one taxa from each defined phylum in the RefSeq-database.

```
sample_ncbi:
    sampling_scheme: "../config/sampling_scheme.ncbi_refseq.yaml"
    database: "RefSeq"
```

__Options for sampling the GTDB Taxonomy__
```
sample_gtdb:
    sampling_scheme: <path>
    completeness: <float>
    contamination: <float>
    gtdb_species_representative: <bool>
    version: <str>
```

`sampling_scheme`: Path to the sampling scheme that will be used. See [Defining a sampling scheme](#defining-a-sampling-scheme) for more details on this.

`completeness`: Exclude taxa with a completeness estimate less than this value (Default: `0`).

`contamination`: Exclude taxa with a contamination estimate larger than this value (Default: `100`).

`gtdb_species_representative`: False will keep all entries while True will only keep entries that are classified as GTDB species representatives.

`version:` Select which version of GTDB to use. E.g. `207` and `214` are supported (Default: `214`).

__Example__

In the example below we will use the sampling scheme defined in `config/sampling_scheme.basic.yaml`. For this example the workflow will sample three taxa for each phylum, but only sample from representative species with an estimated completeness over 90% and an estimated contamination under 5%.

```
sample_gtdb:
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
    taxon: <str>
    version: <str>
```

`bac120`: Number of taxa to sample from the bacterial phylogeny.

`ar53`: Number of taxa to sample from the archaeal phylogeny.

`completeness`: Exclude taxa with a completeness estimate less than this value (Default: `0`).

`contamination`: Exclude taxa with a contamination estimate larger than this value (Default: `100`).

`prune_method`: Select what method to use for pruning, `"shortest"` will keep the taxon with the shortest branch in a leaf-pair, `"longest"` will keep the taxon with the longest branch in a leaf-pair, and `"random"` will randomly select one of the taxa to keep in a leaf-pair (Default: `"shortest"`).

`taxon:` Prune only phylogeny under this taxon, other parts of the phylogeny will be discarded. The taxon must be present in the phylogeny.

`version:` Select which version of GTDB to use, `207` and `214` are supported (Default: `214`).

__Example 1__

In the example below TADA will first remove all taxa with an estimated completeness under 90% and an estimated contamination over 5%. It will then continue to prune the bacterial phylogeny untill 1000 taxa remains. For the archaeal phylogeny it will prune the phylogeny until 200 taxa remains.

```
prune_gtdb:
  bac120: 1000
  ar53: 200
  completeness: 90
  contamination: 5
  prune_method: "shortest"
```

__Example 2__
In the example below TADA will first remove all genomes that are not of high-quality, next it will prune only the Alphaprotobacteria-clade until 100 taxa remains.
```
prune_gtdb:
  bac120: 100
  completeness: 90
  contamination: 5
  prune_method: "shortest"
  taxon: "Alphaprotobacteria"
```


## Defining a sampling scheme

The taxonomic sampling (`sample_gtdb` or `sample_ncbi`) is based on a sampling scheme that is defined in a YAML-file with a structure described below. Examples of sampling schemes can be found in the config-directory.

```
taxonomic_name:
    sampling_level: [domain, phylum, class, order, family, genus, species]
    taxa: <int> or "all"
```
`taxonomic_name`: The name of the taxa to sample from. The name has to be  defined in the GTDB or NCBI taxonomy, depending on which database to sample from. The key-word `all` can be used to sample from all Bacteria and Archaea.

`sampling_level`: The taxonomic level to perform the sampling at.

`taxa`: The number of taxa we want to sample from each group at that taxonomic level. The key-word `all` can also be used, this will keep all taxa in the groups.

Thus, if `taxonomic_name` is set to `Bacteria`, `sampling_level` to `class`, and `taxa` to 3, TADA till sample 3 taxa in each class of the Bacteria.


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
