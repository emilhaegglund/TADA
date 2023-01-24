if config["sample_gtdb"]["gtdb_species_representative"]:
    config["sample_gtdb"]["gtdb_species_representative_opt"] = "--gtdb-representative"
else:
    config["sample_gtdb"]["gtdb_species_representative_opt"] = ""

rule prune_gtdb_phylogeny:
    input:
        phylogeny="gtdb_data/{domain}_r207.tree",
        metadata="gtdb_data/{domain}_metadata_r207.wo_suppressed_records.tsv"
    output:
        phylogeny="prune_gtdb.{domain}.nwk",
        metadata="prune_gtdb.{domain}.metadata.tsv",

    params:
        taxa=lambda wildcards: config["prune_gtdb"]["{}".format(wildcards.domain)],
        completeness=config["prune_gtdb"]["completeness"],
        contamination=config["prune_gtdb"]["contamination"],
    conda:
        "../envs/ete.yaml"
    script:
        "../scripts/prune_gtdb_phylogeny.py"

rule merge_prune_gtdb_output:
    input:
        bacteria_metadata="prune_gtdb.bac120.metadata.tsv",
        archaea_metadata="prune_gtdb.ar53.metadata.tsv"
    output:
        "prune_gtdb.metadata.tsv"
    script:
        """
        ../scripts/merge_tables.py {input.bacteria_metadata} {input.archaea_metadata} {output}
        """

rule subsample_gtdb:
    """
    Use the metadata and taxonomy information to subsample the gtdb-data.
    """
    input:
        metadata="gtdb_data/metadata_r207.wo_suppressed_records.tsv"
    output:
        "sample_gtdb.metadata.tsv"
    params:
        sampling_scheme=config["sample_gtdb"]["sampling_scheme"],
        completeness=config["sample_gtdb"]["completeness"],
        contamination=config["sample_gtdb"]["contamination"],
        gtdb_representative=config["sample_gtdb"]["gtdb_species_representative"]
    script:
        "../scripts/subsample_gtdb_taxonomy.py"

rule download_gtdb_summary:
    input:
        "{method}.metadata.tsv"
    output:
        temp("{method}.ncbi_datasets.tsv")
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        awk -F'\\t'  '{{ print $1 }}' {input} | sed '1d' > accessions.tmp;
        head accessions.tmp;
        datasets summary genome accession --inputfile accessions.tmp \
            --as-json-lines | \
        dataformat tsv genome > {output};
        rm accessions.tmp;
        """

rule filter_gtdb:
    """
    Merge the genome summary files into a single tsv table
    """
    input:
        "{method}.ncbi_datasets.tsv"
    output:
        "{method}.annotation_data.tsv"
    script:
        "../scripts/merge_datasets.py"
