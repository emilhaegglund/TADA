import os

include: "rules/download_gtdb_data.smk"
def all_input(wildcards):
    """
    Generate all inputs for the rule all
    """
    from scripts.define_inputs import gtdb_data_files
    wanted_input = []

    if config["method"]["sample_gtdb"] or config["method"]["trim_gtdb"]:
        wanted_input += gtdb_data_files(config)
    return wanted_input

rule all:
    input: all_input
        #"results/gtdb/bac120_taxonomy_r207.tsv",
        #"results/gtdb/ar53_taxonomy_r207.tsv",
        #"results/gtdb/bac120_metadata_r207.tsv",
        #"results/gtdb/ar53_metadata_r207.tsv",
        #"results/gtdb/sampled_genomes.tsv",
        #directory("results/gtdb/ncbi_proteomes")

### Download required files from GTDB ###

rule sample_gtdb:
    input:
        bac_taxa="results/gtdb/bac120_taxonomy_r207.tsv",
        bac_metadata="results/gtdb/bac120_metadata_r207.tsv"
    output:
        "results/gtdb/sampled_genomes.tsv"
    params:
        taxonomic_level=config["sample_gtdb"]["taxonomic_level"],
        max_taxa=config["sample_gtdb"]["max_taxa"]
    shell:
        """
        python scripts/sample_gtdb.py \
            --gtdb-taxonomy {input.bac_taxa} \
            --gtdb-metadata {input.bac_metadata} \
            --taxonomic-level {params.taxonomic_level} \
            --max-taxa {params.max_taxa} \
            --output {output}
        """
#rule trim_gtdb_taxonomy:
#    input:
#
#    output:
#
#    params:
#        taxa=config["gtdb"]["n_bac_taxa"]
#    conda:
#        "envs/ete.yaml"
#    shell:
#        """
#        python scripts/trim_gtdb_phylogeny.py \
#            --phylogeny {input.phylogeny} \
#            --taxa {params.taxa} \
#            --output {output}
#        """

rule download_proteomes:
    input:
        "results/gtdb/sampled_genomes.tsv"
    output:
        directory("results/gtdb/ncbi_proteomes")
    conda:
        "envs/bit.yaml"
    threads:
        12
    shell:
        """
        mkdir {output};
        cd results/gtdb/;
        input=$(basename {input});
        bit-dl-ncbi-assemblies -w $input -f protein -j {threads};
        mv *.faa.gz {output};
        """

#rule move_not_found:
#    input:
#
#
#
#rule download_genomes:
#    input:
#        "results/gtdb/accessions-not-found.txt"
#    output:
#        directory("results/gtdb/genomes")
#    conda:
#        "envs/bit.yaml"
#    threads:
#       12
#    shell:
#        """
#        bit-dl-ncbi-assemblies -w {input} -f fasta -j {threads}
#        """
#
##TODO: What to do with suppresed genomes?
#
rule run_prodigal:
    input:
        "results/gtdb/"
    output:
        "results/gtdb/prodigal_annotation"
    conda:
        "envs/prodigal.yaml"
    shell:
        "prodigal"

#rule download_gtdb_format_taxonomy:
#
#rule create_prot2taxa_map:
#
#rule build_diamond:
#    input:
#
#    output:
#
#    conda:
#        "envs/diamond.yaml"
