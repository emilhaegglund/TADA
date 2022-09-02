rule sample_gtdb:
    """
    Use the metadata and taxonomy information to subsample the gtdb-data.
    """
    input:
        bac_taxa=config["paths"]["results"] + "/gtdb_data/{domain}_taxonomy_r207.tsv",
        bac_metadata=config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.tsv"
    output:
        genbank=config["paths"]["results"] + "/gtdb_sampling/sampled_{domain}_genbank_genomes.tsv",
        refseq=config["paths"]["results"] + "/gtdb_sampling/sampled_{domain}_refseq_genomes.tsv"
    params:
        taxonomic_level=config["sample_gtdb"]["taxonomic_level"],
        max_taxa=config["sample_gtdb"]["max_taxa"],
        completeness=config["sample_gtdb"]["completeness"],
        contamination=config["sample_gtdb"]["contamination"],
        gtdb_representative=config["sample_gtdb"]["gtdb_representative"]
    shell:
        """
        python scripts/sample_gtdb.py \
            --gtdb-taxonomy {input.bac_taxa} \
            --gtdb-metadata {input.bac_metadata} \
            --taxonomic-level {params.taxonomic_level} \
            --max-taxa {params.max_taxa} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --gtdb-representative {params.gtdb_representative} \
            --output-refseq {output.refseq} \
            --output-genbank {output.genbank}
        """

rule ncbi_download_genomes_refseq:
    input:
        config["paths"]["results"] + "/gtdb_sampling/sampled_{domain}_{database}_genomes.tsv"
    output:
        dir(config["paths"]["results"] + "/gtdb_sampling/sampled_{domain}_{database}_genomes/")
    params:
        db='{database}'
    threads:
        12
    conda:
        "../envs/ncbi-download-genomes.yaml"
    shell:
        """
        ncbi-download-genomes -F protein-fasta -A {input} --flat-output -p {threads} -o {output} -s {params.db} all
        """

rule trim_gtdb_taxonomy:
    input:
        phylogeny=config["paths"]["results"] + "/gtdb_data/bac120_r207.tree"
    output:
        config["paths"]["results"] + "/gtdb_trimming/bac120_r207_trimmed_genomes.tsv"
    params:
        taxa=config["trim_gtdb"]["number_bac"]
    conda:
        "../envs/ete.yaml"
    shell:
        """
        python scripts/trim_gtdb_phylogeny.py \
            --phylogeny {input.phylogeny} \
            --taxa {params.taxa} \
            --output {output}
        """
