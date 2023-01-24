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
    return expand(os.path.join("ncbi_proteomes",
                               "{accession}.faa"),
                               accession=accession
                               )

checkpoint get_genomes_w_annotation:
    input:
        config["method"] + ".annotation_data.tsv"
    output:
        directory("accession_with_annotation")
    params:
        status="with-annotation"
    script:
        "../scripts/divide_accessions_on_annotation_status.py"

rule download_proteomes:
    input:
        get_genomes_with_annotation
    output:
        temp("ncbi_proteomes/{accession}.zip")
    params:
        acc="{accession}"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession {params.acc}\
            --include protein \
            --filename {output} \
            --no-progressbar
        """

rule download_genomes:
    input:
        get_genomes_with_annotation
    output:
        temp("ncbi_genomes/{accession}.zip")
    params:
        acc="{accession}"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession {params.acc}\
            --include genome \
            --filename {output} \
            --no-progressbar
        """

rule unzip_ncbi_proteomes:
    input:
        "ncbi_proteomes/{accession}.zip"
    output:
        "ncbi_proteomes/{accession}.faa"
    params:
        acc="{accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/*.faa > {output}"

rule unzip_ncbi_genomes_with_annotation:
    input:
        "ncbi_genomes/{accession}.zip"
    output:
        "ncbi_genomes/{accession}.faa"
    params:
        acc="{accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/*.faa > {output}"

# Genomes without annotation
def get_genomes_without_annotation(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join(ck_output,
            "{prokka_accession}.txt"),
            prokka_accession=accession
                               )
def get_prokka_proteomes(wildcards):
    ck_output = checkpoints.get_genomes_wo_annotation.get(**wildcards).output[0]
    accession, = glob_wildcards(os.path.join(ck_output, "{sample}.txt"))
    return expand(os.path.join("prokka",
                               "{prokka_accession}",
                               "{prokka_accession}.faa"),
                               prokka_accession=accession
                               )

checkpoint get_genomes_wo_annotation:
    input:
        config["method"] + ".annotation_data.tsv"
    output:
        directory("accession_without_annotation")
    params:
        status="without-annotation"
    script:
        "../scripts/divide_accessions_on_annotation_status.py"

rule download_genomes_wo_annotation:
    input:
        get_genomes_without_annotation
    output:
        temp("genomes_wo_annotation/{prokka_accession}.zip")
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
        "genomes_wo_annotation/{prokka_accession}.zip"
    output:
        "genomes_wo_annotation/{prokka_accession}.fna"
    params:
        acc="{prokka_accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/*.fna > {output}"

rule prokka:
    input:
        "genomes_wo_annotation/{prokka_accession}.fna"
    output:
        "prokka/{prokka_accession}/{prokka_accession}.faa"
    params:
        prefix="{prokka_accession}",
        outdir="prokka/{prokka_accession}"
    threads:
        4
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir} --prefix {params.prefix} --cpus {threads} {input} --force
        """

rule collect_proteomes:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        directory("proteomes/")
    shell:
        """
        mkdir {output}
        for prot in {input.prokka} {input.ncbi}
        do
            ln -s $prot {output}
        done
        """

rule build_diamond:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        "diamond_db/sdbw.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    shell:
        """
        cat {input.prokka} {input.ncbi} | diamond makedb --db {output} -p {threads}
        """

rule make_blast_db:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        "ncbi_blastp_db/sdbw.pdb"
    params:
        prefix="ncbi_blastp_db/sdbw",
        title="sdbw"
    conda:
        "../envs/ncbi_blast.yaml"
    shell:
        """
        cat {input.prokka} {input.ncbi} | makeblastdb -in - -dbtype prot -out {params.prefix} -title {params.title};
        """
