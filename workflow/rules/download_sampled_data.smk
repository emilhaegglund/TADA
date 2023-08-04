config["include"] = []
if config["downloads"]["proteomes"]:
    config["include"].append("protein")
if config["downloads"]["genomes"]:
    config["include"].append("genome")
if config["downloads"]["cds"]:
    config["include"].append("cds")
config["include"] = ",".join(config["include"])

# Genomes wiht annotation
def get_genomes_with_annotation(wildcards):
    ck_output = checkpoints.get_genomes_w_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join(ck_output,
            "{accession}.txt"),
            accession=accession
                               )

def get_ncbi_proteomes(wildcards):
    ck_output = checkpoints.get_genomes_w_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("workflow_files",
                               "ncbi_proteomes",
                               "{accession}.faa"),
                               accession=accession
                               )

def get_ncbi_cds(wildcards):
    ck_output = checkpoints.get_genomes_w_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("workflow_files",
                               "ncbi_cds",
                               "{accession}.ffn"),
                               accession=accession
                               )

def get_ncbi_genomes(wildcards):
    ck_output = checkpoints.get_genomes_w_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("workflow_files",
                               "ncbi_genomes",
                               "{accession}.fna"),
                               accession=accession
                               )

checkpoint get_genomes_w_annotation:
    input:
        config["method"] + ".annotation_data.tsv"
    output:
        directory("workflow_files/accession_with_annotation")
    params:
        status="with-annotation"
    script:
        "../scripts/divide_accessions_on_annotation_status.py"

rule download_annotated_data:
    output:
        temp("workflow_files/ncbi_annotated_data/{accession}.zip")
    params:
        acc="{accession}",
        inc=config["include"]
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession {params.acc}\
            --include {params.inc} \
            --filename {output} \
            --no-progressbar
        """

rule unzip_ncbi_proteomes:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_proteomes/{accession}.faa"
    params:
        acc="{accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/protein.faa > {output}"

rule unzip_ncbi_cds:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_cds/{accession}.ffn"
    params:
        acc="{accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/cds_from_genomic.fna > {output}"

rule unzip_ncbi_genomes_with_annotation:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_genomes/{accession}.fna"
    params:
        acc="{accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/{params.acc}*.fna > {output}"

# Genomes without annotation
def get_genomes_without_annotation(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join(ck_output,
            "{prokka_accession}.txt"),
            prokka_accession=accession
                               )

def get_genomes_without_annotation_for_genome_dir(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("workflow_files", "genomes_wo_annotation",
            "{accession}.fna"),
            accession=accession
                               )
def get_prokka_proteomes(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("workflow_files",
                               "prokka",
                               "{prokka_accession}",
                               "{prokka_accession}.faa"),
                               prokka_accession=accession
                               )
def get_prokka_cds(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("workflow_files",
                               "prokka",
                               "{prokka_accession}",
                               "{prokka_accession}.ffn"),
                               prokka_accession=accession
                               )


checkpoint get_genomes_wo_annotation:
    input:
        config["method"] + ".annotation_data.tsv"
    output:
        directory("workflow_files/accession_without_annotation")
    params:
        status="without-annotation"
    script:
        "../scripts/divide_accessions_on_annotation_status.py"

rule download_genomes_wo_annotation:
    output:
        temp("workflow_files/genomes_wo_annotation/{prokka_accession}.zip")
    params:
        acc="{prokka_accession}"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession {params.acc}\
            --include genome \
            --filename {output} \
            --no-progressbar
        """

rule unzip_genomes_wo_annotation:
    input:
        "workflow_files/genomes_wo_annotation/{prokka_accession}.zip"
    output:
        "workflow_files/genomes_wo_annotation/{prokka_accession}.fna"
    params:
        acc="{prokka_accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/*.fna > {output}"

rule prokka:
    input:
        "workflow_files/genomes_wo_annotation/{prokka_accession}.fna"
    output:
        proteomes="workflow_files/prokka/{prokka_accession}/{prokka_accession}.faa",
        cds="workflow_files/prokka/{prokka_accession}/{prokka_accession}.ffn"
    params:
        prefix="{prokka_accession}",
        outdir="workflow_files/prokka/{prokka_accession}"
    threads:
        4
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir} --prefix {params.prefix} --cpus {threads} {input} --force
        """

# Use a script to link files from the NCBI and Prokka directories to
# a common directory. The script is used to get the absolute paths
# to the linked files for better integration with downstream analysis tools.
rule collect_proteomes:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        directory("proteomes/")
    script:
        "../scripts/link_files.py"

rule collect_cds:
    input:
        prokka=get_prokka_cds,
        ncbi=get_ncbi_cds,
    output:
        directory("cds/")
    script:
        "../scripts/link_files.py"

rule collect_genomes:
    input:
        ncbi=get_ncbi_genomes,
        also_ncbi=get_genomes_without_annotation_for_genome_dir
    output:
        directory("genomes/")
    script:
        "../scripts/link_files.py"


rule build_diamond_protein:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        "databases/diamond_proteomes_db/sdbw_proteomes.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    shell:
        """
        cat {input.prokka} {input.ncbi} | diamond makedb --db {output} -p {threads}
        """

rule build_blast_protein:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        "databases/blast_proteomes_db/sdbw_proteomes.pdb"
    params:
        prefix="databases/blast_proteomes_db/sdbw_proteomes",
        title="sdbw_proteomes"
    conda:
        "../envs/ncbi_blast.yaml"
    shell:
        """
        cat {input.prokka} {input.ncbi} | makeblastdb -in - -dbtype prot -out {params.prefix} -title {params.title};
        """

rule build_blast_cds:
    input:
        prokka=get_prokka_cds,
        ncbi=get_ncbi_cds
    output:
        "databases/blast_cds_db/sdbw_cds.ndb"
    params:
        prefix="databases/blast_cds_db/sdbw_cds",
        title="sdbw_cds"
    conda:
        "../envs/ncbi_blast.yaml"
    shell:
        """
        cat {input.prokka} {input.ncbi} | makeblastdb -in - -dbtype nucl -out {params.prefix} -title {params.title}
        """

rule build_blast_genome:
    input:
        ncbi=get_ncbi_genomes,
        also_ncbi=get_genomes_without_annotation_for_genome_dir
    output:
        "databases/blast_genomes_db/sdbw_genomes.ndb"
    params:
        prefix="databases/blast_genomes_db/sdbw_genomes",
        title="sdbw_genomes"
    conda:
        "../envs/ncbi_blast.yaml"
    shell:
        """
        cat {input.ncbi} {input.also_ncbi} | makeblastdb -in - -dbtype nucl -out {params.prefix} -title {params.title}
        """
