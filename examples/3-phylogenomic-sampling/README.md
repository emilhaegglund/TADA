In this example TADA has been configured to use the phylogenomic sampling approach to sample 100 MAGs/genomes from the Alphaproteobacteria class by pruning the bacterial phylogeny provided by GTDB r214.

This example can be run using the following command:
```
snakemake -s ../../workflow/Snakefile --use-conda --conda-frontend mamba -j4 --configfile phylogenomic-config.yaml
```
