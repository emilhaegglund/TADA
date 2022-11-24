if config["sample_gtdb"]["gtdb_species_representative"]:
    config["sample_gtdb"]["gtdb_species_representative_opt"] = "--gtdb-representative"
else:
    config["sample_gtdb"]["gtdb_species_representative_opt"] = ""

rule prune_gtdb_phylogeny:
    input:
        phylogeny=config["base_dir"] + "/gtdb_data/{domain}_r207.tree",
        metadata=config["base_dir"] + "/gtdb_data/{domain}_metadata_r207.wo_suppressed_records.tsv"
    output:
        metadata=config["base_dir"] + "/{prefix}.prune_gtdb.{domain}.metadata.tsv",
        tree=config["base_dir"] + "/{prefix}.prune_gtdb.{domain}.nwk"

    params:
        taxa=lambda wildcards: config["prune_gtdb"]["{}".format(wildcards.domain)],
        completeness=config["prune_gtdb"]["completeness"],
        contamination=config["prune_gtdb"]["contamination"],
    conda:
        "../envs/ete.yaml"
    shell:
        """
        python scripts/prune_gtdb_phylogeny.py \
            --phylogeny {input.phylogeny} \
            --gtdb-metadata {input.metadata} \
            --taxa {params.taxa} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --output-metadata {output.metadata} \
            --output-tree {output.tree} \
        """

rule merge_prune_gtdb_output:
    input:
        bacteria_metadata=config["base_dir"] + "/{prefix}.prune_gtdb.bac120.metadata.tsv",
        archaea_metadata=config["base_dir"] + "/{prefix}.prune_gtdb.ar53.metadata.tsv"
    output:
        config["base_dir"] + "/{prefix}.prune_gtdb.metadata.tsv"
    shell:
        """
        python scripts/merge_tables.py {input.bacteria_metadata} {input.archaea_metadata} {output}
        """

rule subsample_gtdb:
    """
    Use the metadata and taxonomy information to subsample the gtdb-data.
    """
    input:
        metadata=config["base_dir"] + "/gtdb_data/metadata_r207.wo_suppressed_records.tsv"
    output:
        config["base_dir"] + "/{prefix}.sample_gtdb.metadata.tsv"
    params:
        sampling_scheme=config["sample_gtdb"]["sampling_scheme"],
        completeness=config["sample_gtdb"]["completeness"],
        contamination=config["sample_gtdb"]["contamination"],
        gtdb_representative=config["sample_gtdb"]["gtdb_species_representative_opt"]
    shell:
        """
        python scripts/subsample_gtdb.py \
            --gtdb-metadata {input.metadata} \
            --sampling-scheme {params.sampling_scheme} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            {params.gtdb_representative} \
            --output {output}
        """
