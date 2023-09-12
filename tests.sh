#!/usr/bin/env bash
# Run through a set of different configuration-files to test the workflow
set -e
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.sample_refseq.yaml -j 12;
echo "Sample RefSeq finished"

snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.sample_genbank.yaml -j 12;
echo "Sample GenBank finished"

snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.sample_gtdb.yaml -j 12;
echo "Sample GTDB finished";

snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.prune_gtdb.yaml -j 12;
echo "Prune GTDB finished";

snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.prune_gtdb_taxon.yaml -j 12;
echo "Prune GTDB with taxon option finished";

snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.sample_gtdb_required.yaml -j 12;
echo "Sample GTDB with required genomes finished";

snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --conda-prefix ../conda-envs --configfile config/config.sample_refseq_required.yaml -j 12;
echo "Sample NCBI with required genomes finsihed";
