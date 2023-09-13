In this example TADA have been configured to sample 1 high-quality MAG/genome for each family in the Alphaproteobacteria order. This example will not download any genomic data, only perform the sampling.

It can be run using the following command
```
snakemake -s ../../workflow/Snakefile --use-conda --conda-frontend mamba -j4 --configfile basic-config.yaml
```
