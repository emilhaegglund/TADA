
if "prune_gtdb" not in config.keys():
    config["prune_gtdb"] = {"bac120": 0,
                            "ar53": 0,
                            "completeness": 0,
                            "contamination": 100
                            }

if "sample_gtdb" not in config.keys():
    config["sample_gtdb"] = {"sampling_scheme": "",
                            "completeness": 0,
                            "contamination": 100,
                            "gtdb_species_representative": False
                            }

validate(config, schema="../validation_schemes/config.schema.yaml")

if config["sample_gtdb"]["gtdb_species_representative"]:
    config["sample_gtdb"]["gtdb_species_representative_opt"] = "--gtdb-representative"
else:
    config["sample_gtdb"]["gtdb_species_representative_opt"] = ""

rule prune_gtdb_phylogeny:
    input:
        phylogeny="gtdb_data/{domain}_r{version}.tree",
        metadata="gtdb_data/metadata_r{version}.wo_suppressed_records.tsv"
    output:
        phylogeny="prune_gtdb.{domain}.nwk",
        metadata="prune_gtdb.{domain}.metadata.tsv",
    params:
        taxa=lambda wildcards: config["prune_gtdb"]["{}".format(wildcards.domain)],
        completeness=config["prune_gtdb"]["completeness"],
        contamination=config["prune_gtdb"]["contamination"],
        seed=config["seed"]
    conda:
        "../envs/ete.yaml"
    script:
        "../scripts/prune_gtdb_phylogeny.py"

rule merge_prune_gtdb_output:
    input:
        bacteria_metadata="prune_gtdb.bac120.metadata.tsv",
        archaea_metadata="prune_gtdb.ar53.metadata.tsv"
    output:
        metadata="prune_gtdb.metadata.tsv"
    script:
        "../scripts/merge_pruned_tables.py"

def subsample_input(wildcards):
    return "gtdb_data/metadata_r" + str(config["sample_gtdb"]["version"]) + ".wo_suppressed_records.tsv"

rule subsample_gtdb:
    """
    Use the metadata and taxonomy information to subsample the gtdb-data.
    """
    input:
        metadata=subsample_input
        #metadata="gtdb_data/metadata_r{version}.wo_suppressed_records.tsv"
    output:
        "sample_gtdb.metadata.tsv"
    params:
        sampling_scheme=config["sample_gtdb"]["sampling_scheme"],
        completeness=config["sample_gtdb"]["completeness"],
        contamination=config["sample_gtdb"]["contamination"],
        gtdb_representative=config["sample_gtdb"]["gtdb_species_representative"],
        seed=config["seed"]
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
