In this example TADA is configured to perform hierarchical taxonomic sampling of the Rhizobiales order in GTDB r214. It will sample 10 MAGs/genomes from Bartonella bacilliformis, 2 MAGs/genomes from each species in the Bartonella genus, and 5 MAGs/genomes from the Rhizobiales order, except from the Bartonella genus. This example will also download the proteomes for these MAGs/genomes and build a Diamond-database.

To run the example, use the command below. This will use 4 cores, to speed up things this number can be increased if you have access to more cores.
```
snakemake -s ../../workflow/Snakefile --use-conda --conda-frontend mamba -j4 --configfile advanced-config.yaml
```
