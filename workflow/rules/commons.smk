import yaml
import os


def get_database_path(config, db, db_type, directory, file_name, wanted_input):
    if config["databases"].get(db):
        if config["downloads"].get(db_type):
            wanted_input.append(os.path.join("databases", directory, file_name))
        else:
            print(f"Building {db} database requires turning the download {db_type} to True")

def all_input(wildcards):
    """
    Generate all inputs for the rule all
    """
    wanted_input = []
    base_dir = config["workdir"]

    if config["method"] == "sample_ncbi":
        with open(config["sample_ncbi"]["sampling_scheme"], "r") as stream:
            sampling_scheme = yaml.safe_load(stream)
        for taxa in sampling_scheme.keys():
            if taxa == 'all':  # If all, only include Bacteria and Archaea
                print('Bact, Arch')
                wanted_input.append(os.path.join("ncbi_data", "taxa", "Bacteria.tsv"))
                wanted_input.append(os.path.join("ncbi_data", "taxa", "Archaea.tsv"))
            else:
                print(taxa)
                wanted_input.append(os.path.join("ncbi_data", "taxa", taxa.replace(" ", "_") + ".tsv"))
        wanted_input.append(os.path.join("ncbi_data", "datasets.tsv"))
        wanted_input.append("sample_ncbi.annotation_data.tsv")

    if config["method"] == "sample_gtdb":
        #wanted_input.append("sample_gtdb.annotation_data.tsv")
        wanted_input.append("sample_gtdb.metadata.tsv")

    if config["method"] == "prune_gtdb":
        gtdb_version = str(config["prune_gtdb"]["version"])
        wanted_input.append("prune_gtdb.metadata.tsv")

    # Handle download options
    if "downloads" in config:
        if "genomes" in config["downloads"]:
            if config["downloads"]["genomes"]:
                wanted_input.append("genomes")
        if "cds" in config["downloads"]:
            if config["downloads"]["cds"]:
                wanted_input.append("cds")
        if "proteomes" in config["downloads"]:
            if config["downloads"]["proteomes"]:
                wanted_input.append("proteomes")
        if "gff3" in config["downloads"]:
            if config["downloads"]["gff3"]:
                wanted_input.append("gff3")

    # Handle blast databases options
    if "databases" in config:
        get_database_path(config, "diamond_protein", "proteomes", "diamond_proteomes_db", "tada_proteomes.dmnd", wanted_input)
        get_database_path(config, "blast_protein", "proteomes", "blast_proteomes_db", "tada_proteomes.pdb", wanted_input)
        get_database_path(config, "blast_genomes", "genomes", "blast_genomes_db", "tada_genomes.ndb", wanted_input)
        get_database_path(config, "blast_cds", "cds", "blast_cds_db", "tada_cds.ndb", wanted_input)

    return wanted_input

def get_accession_w_annotation(wildcards):
    ck_output = checkpoints.get_genomes_w_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return(accession)

def get_accession_wo_annotation(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return(accession)

def get_genomes_with_annotation(wildcards):
    ck_output = checkpoints.get_genomes_w_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join(ck_output, "{accession}.txt"), accession=accession)

def get_ncbi_proteomes(wildcards):
    accession = get_accession_w_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "ncbi_proteomes",
                               "{accession}.faa"),
                               accession=accession
                               )

def get_ncbi_cds(wildcards):
    accession = get_accession_w_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "ncbi_cds",
                               "{accession}.ffn"),
                               accession=accession
                               )

def get_ncbi_genomes(wildcards):
    accession = get_accession_w_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "ncbi_genomes",
                               "{accession}.fna"),
                               accession=accession
                               )

def get_ncbi_gff3(wildcards):
    accession = get_accession_w_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "ncbi_gff3",
                               "{accession}.gff"),
                               accession=accession
                               )
# Genomes without annotation
def get_genomes_without_annotation(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join(ck_output,
            "{prokka_accession}.txt"),
            prokka_accession=accession
                               )

def get_genomes_without_annotation_for_genome_dir(wildcards):
    accession = get_accession_wo_annotation(wildcards)
    return expand(os.path.join("workflow_files", "genomes_wo_annotation",
            "{accession}.fna"),
            accession=accession
                               )
def get_prokka_proteomes(wildcards):
    accession = get_accession_wo_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "prokka",
                               "{prokka_accession}",
                               "{prokka_accession}.faa"),
                               prokka_accession=accession
                               )
def get_prokka_cds(wildcards):
    accession = get_accession_wo_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "prokka",
                               "{prokka_accession}",
                               "{prokka_accession}.ffn"),
                               prokka_accession=accession
                               )

def get_prokka_gff3(wildcards):
    accession = get_accession_wo_annotation(wildcards)
    return expand(os.path.join("workflow_files",
                               "prokka",
                               "{prokka_accession}",
                               "{prokka_accession}.gff"),
                               prokka_accession=accession
                               )
# gtdb_taxonomy_phylogeny
def bac120_prune_input(wildcards):
    return "prune_gtdb.bac120_r" + str(config["prune_gtdb"]["version"]) + ".metadata.tsv"

def ar53_prune_input(wildcards):
    return "prune_gtdb.ar53_r" + str(config["prune_gtdb"]["version"]) + ".metadata.tsv"

def subsample_input(wildcards):
    return "gtdb_data/metadata_r" + str(config["sample_gtdb"]["version"]) + ".wo_suppressed_records.tsv"
