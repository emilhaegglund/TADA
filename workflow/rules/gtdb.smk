rule prune_gtdb_phylogeny:
    input:
        phylogeny=config["paths"]["results"] + "/gtdb_data/{domain}_r207.tree",
        metadata=config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.wo_suppressed_records.tsv"
    output:
        refseq=config["paths"]["results"] + "/prune_gtdb/sampled_{domain}_refseq_accessions.tsv",
        genbank=config["paths"]["results"] + "/prune_gtdb/sampled_{domain}_genbank_accessions.tsv",
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
            --taxa {params.taxa} \
            --gtdb-metadata {input.metadata} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --output-refseq {output.refseq} \
            --output-genbank {output.genbank}
        """

rule subsample_gtdb:
    """
    Use the metadata and taxonomy information to subsample the gtdb-data.
    """
    input:
        metadata=config["paths"]["results"] + "/gtdb_data/metadata_r207.wo_suppressed_records.tsv"
    output:
        config["paths"]["results"] + "/subsample_gtdb/sampled_accessions.metadata.tsv",
    params:
        sampling_scheme=config["subsample_gtdb"]["sampling_scheme"],
        completeness=config["subsample_gtdb"]["completeness"],
        contamination=config["subsample_gtdb"]["contamination"],
        gtdb_representative=config["subsample_gtdb"]["gtdb_species_representative"]
    shell:
        """
        python scripts/subsample_gtdb.py \
            --gtdb-metadata {input.metadata} \
            --sampling-scheme {params.sampling_scheme} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --gtdb-representative {params.gtdb_representative} \
            --output {output}
        """

rule get_sampled_accessions:
    input:
        config["paths"]["results"] + "/subsample_gtdb/sampled_accessions.metadata.tsv",
    output:
        config["paths"]["results"] + "/subsample_gtdb/sampled_accessions.txt",
    shell:
        "cut -f1 {input} | grep -v 'Accessions' > {output}"
