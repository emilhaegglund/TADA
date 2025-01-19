config["include"] = []
if "downloads" in config:
    if "proteomes" in config["downloads"]:
        config["include"].append("protein")
    if "genomes" in config["downloads"]:
        config["include"].append("genome")
    if "cds" in config["downloads"]:
        config["include"].append("cds")
    if "gff3" in config["downloads"]:
        config["include"].append("gff3")
    config["include"] = ",".join(config["include"])


checkpoint get_genomes_w_annotation:
    input:
        config["method"] + ".annotation_data.tsv"
    output:
        directory("workflow_files/accession_with_annotation")
    params:
        status="with-annotation"
    conda:
        "../envs/base.yaml"
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
    log:
        "logs/download_annotated_data.{accession}.log"
    shell:
        """
        datasets download genome accession {params.acc}\
            --include {params.inc} \
            --filename {output} \
            --no-progressbar 2> {log}
        """

rule unzip_ncbi_proteomes:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_proteomes/{accession}.faa"
    params:
        acc="{accession}"
    conda:
        "../envs/base.yaml"
    log:
        "logs/unzip_ncbi_proteomes.{accession}.log"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/protein.faa > {output}"

rule unzip_ncbi_cds:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_cds/{accession}.ffn"
    params:
        acc="{accession}"
    conda:
        "../envs/base.yaml"
    log:
        "logs/unzip_ncbi_cds.{accession}.log"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/cds_from_genomic.fna > {output} 2> {log}"

rule unzip_ncbi_gff3:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_gff3/{accession}.gff"
    params:
        acc="{accession}"
    conda:
        "../envs/base.yaml"
    log:
        "logs/unzip_ncbi_gff3.{accession}.log"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/genomic.gff > {output} 2> {log}"

rule unzip_ncbi_genomes_with_annotation:
    input:
        "workflow_files/ncbi_annotated_data/{accession}.zip"
    output:
        "workflow_files/ncbi_genomes/{accession}.fna"
    params:
        acc="{accession}"
    conda:
        "../envs/base.yaml"
    log:
        "logs/unzip_ncbi_genomes_with_annotation.{accession}.log"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/{params.acc}*.fna > {output} 2> {log}"

checkpoint get_genomes_wo_annotation:
    input:
        config["method"] + ".annotation_data.tsv"
    output:
        directory("workflow_files/accession_without_annotation")
    conda:
        "../envs/base.yaml"
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
    log:
        "logs/download_genomes_wo_annotation.{prokka_accession}.log"
    shell:
        """
        datasets download genome accession {params.acc}\
            --include genome \
            --filename {output} \
            --no-progressbar 2> {log}
        """

rule unzip_genomes_wo_annotation:
    input:
        "workflow_files/genomes_wo_annotation/{prokka_accession}.zip"
    output:
        "workflow_files/genomes_wo_annotation/{prokka_accession}.fna"
    conda:
        "../envs/biopython.yaml"
    params:
        acc="{prokka_accession}"
    log:
        "logs/unzip_genomes_wo_annotation.{prokka_accession}.log"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/*.fna > {output} 2> {log}"

rule prokka:
    input: 
        "workflow_files/genomes_wo_annotation/{prokka_accession}.fna"
    output:
        proteomes="workflow_files/prokka/{prokka_accession}/{prokka_accession}.faa",
        cds="workflow_files/prokka/{prokka_accession}/{prokka_accession}.ffn",
        gff3="workflow_files/prokka/{prokka_accession}/{prokka_accession}.gff"
    params:
        prefix="{prokka_accession}",
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads:
        4
    conda: 
        "../envs/prokka.yaml"
    log:
        "logs/prokka_{prokka_accession}.log"
    shell:
        """
        prokka --outdir {params.outdir} --force --prefix {params.prefix} --cpus {threads} {input} 2> {log}
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
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/link_files.py"

rule collect_cds:
    input:
        prokka=get_prokka_cds,
        ncbi=get_ncbi_cds,
    output:
        directory("cds/")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/link_files.py" 

rule collect_genomes:
    input:
        ncbi=get_ncbi_genomes,
        also_ncbi=get_genomes_without_annotation_for_genome_dir
    output:
        directory("genomes/")
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/link_files.py"

rule collect_gff3:
    input:
        prokka=get_prokka_gff3,
        ncbi=get_ncbi_gff3,
    output:
        directory("gff3/")
    conda: "../envs/base.yaml"
    script:
        "../scripts/link_files.py"

rule build_diamond_protein:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        "databases/diamond_proteomes_db/tada_proteomes.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    log:
      "logs/build_diamond_proteins.log"
    shell:
        """
        cat {input.prokka} {input.ncbi} | diamond makedb --db {output} -p {threads} 2> {log}
        """

rule build_blast_protein:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes,
    output:
        "databases/blast_proteomes_db/tada_proteomes.pdb"
    params:
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        title="tada_proteomes"
    conda:
        "../envs/ncbi-blast.yaml"
    log:
        "logs/build_blast_protein.log"
    shell:
        """
        cat {input.prokka} {input.ncbi} | \
        seqkit rmdup | \
        makeblastdb -in - -dbtype prot -out {params.prefix} -title {params.title} -parse_seqids; 2> {log}
        """

rule build_blast_cds:
    input:
        prokka=get_prokka_cds,
        ncbi=get_ncbi_cds
    output:
        "databases/blast_cds_db/tada_cds.ndb"
    params:
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        title="tada_cds"
    conda:
        "../envs/ncbi-blast.yaml"
    log:
        "logs/build_blast_cds.log"
    shell:
        """
        cat {input.prokka} {input.ncbi} | \
        seqkit rmdup | \
        makeblastdb -in - -dbtype nucl -out {params.prefix} -title {params.title} -parse_seqids; 2> {log}
        """

rule build_blast_genome:
    input:
        ncbi=get_ncbi_genomes,
        also_ncbi=get_genomes_without_annotation_for_genome_dir
    output:
        "databases/blast_genomes_db/tada_genomes.ndb"
    params:
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        title="tada_genomes"
    conda:
        "../envs/ncbi-blast.yaml"
    log:
        "logs/build_blast_genome.log"
    shell:
        """
        cat {input.ncbi} {input.also_ncbi} | \
        seqkit rmdup | \
        makeblastdb -in - -dbtype nucl -out {params.prefix} -title {params.title} -parse_seqids; 2> {log}
        """