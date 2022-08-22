rule sample_gtdb:
    """
    Use the metadata and taxonomy information to subsample the gtdb-data.
    """
    input:
        bac_taxa=config["paths"]["results"] + "/gtdb_data/{domain}_taxonomy_r207.tsv",
        bac_metadata=config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.tsv"
    output:
        config["paths"]["results"] + "/gtdb_sampling/sampled_{domain}_genomes.tsv"
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
            --output {output}
        """

rule trim_gtdb_taxonomy:
    input:
        phylogeny=config["paths"]["results"] + "/gtdb_data/bac120_r207.tree"
    output:
        config["paths"]["results"] + "/gtdb_trimming/bac120_r207_trimmed_genomes.tsv"
    params:
        taxa=config["trim_gtdb"]["number_bac"]
    conda:
        "envs/ete.yaml"
    shell:
        """
        python scripts/trim_gtdb_phylogeny.py \
            --phylogeny {input.phylogeny} \
            --taxa {params.taxa} \
            --output {output}
        """
