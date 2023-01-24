def get_prokka_proteomes(wildcards):
    ck_output = checkpoints.get_genomes_wo_annoation.get(**wildcards).output[0]
    prokka_accession, = glob_wildcards(os.path.join(ck_output, "{accession}.txt"))
    return expand(os.path.join(config["base_dir"], "prokka", "{prokka_accession}", "{prokka_accession}.faa"), prokka_accession=prokka_accession)

def get_ncbi_proteomes(wildcards):
    ncbi_proteomes = checkpoints.unzip_proteomes.get(**wildcards).output[0]
    accession = glob_wildcards(ncbi_proteomes + "/ncbi_dataset/data/{accession}/proteins.faa")
    return expand(os.path.join(config["base_dir"], "ncbi_dataset", "data", "{accession}", "proteins.faa"), accession=accession)

rule get_sampled_accessions:
    input:
       config["base_dir"] + "/{prefix}.{method}.metadata.tsv"
    output:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt",
    shell:
        "cut -f1 {input} | grep -v 'assembly_accession' > {output}"

rule get_refseq_accessions:
    input:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt",
    output:
        temporary(config["base_dir"] + "/{prefix}.{method}.refseq.accessions.txt")
    shell:
        """
        grep 'GCF' {input} > {output} || true;
        """

rule get_genbank_accessions:
    input:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt"
    output:
        temporary(config["base_dir"] + "/{prefix}.{method}.genbank.accessions.txt")
    shell:
        """
        echo 'accession' > {output};
        grep 'GCA' {input} >> {output} || true;
        """

rule datasets_download_proteome:
    input:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt"
    output:
        temp(config["base_dir"] + "/{prefix}.{method}.proteome_data.zip")
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession --inputfile {input} \
            --include protein \
            --filename {output} \
            --no-progressbar
        """

rule datasets_download_genome:
    input:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt"
    output:
        temp(config["base_dir"] + "/{prefix}.{method}.genome_data.zip")
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        datasets download genome accession --inputfile {input} \
            --include genome \
            --filename {output} \
            --no-progressbar
        """

checkpoint unzip_proteomes:
    input:
        config["base_dir"] + "/{prefix}.{method}.proteome_data.zip"
    output:
        directory(config["base_dir"] + "/{prefix}.{method}.proteome_data/")
    shell:
        "unzip -d {output} {input}"

rule unzip_genomes:
    input:
        config["base_dir"] + "/{prefix}.{method}.genome_data.zip"
    output:
        directory(config["base_dir"] + "/{prefix}.{method}.genome_data")
    shell:
        "unzip -d {output} {input}"

checkpoint get_genomes_wo_annoation:
    input:
        config["base_dir"] + "/{prefix}.{method}.metadata.tsv"
    output:
        directory(config["base_dir"] + "/{prefix}.{method}.accessions_wo_annotation/")
    shell:
        """
        mkdir {output};
        python scripts/extract_genome_wo_annotation.py {input} {output};
        """

rule unzip_genome_wo_annotation:
    input:
        config["base_dir"] + "/{prefix}.{method}.genomes_wo_annotation_zipped/{accession}.zip"
    output:
        config["base_dir"] + "/{prefix}.{method}.genomes_wo_annotation/{accession}.fna"
    params:
        acc="{accession}"
    shell:
        "unzip -p {input} ncbi_dataset/data/{params.acc}/*.fna > {output}"

rule prokka:
    input:
        fasta=config["base_dir"] + "{prefix}.{method}.genomes_wo_annotation/{prokka_accession}.fna",
    output:
        config["base_dir"] + "/prokka/{prokka_accession}/{prokka_accession}.faa"
    params:
        prefix="{prokka_accession}",
        outdir=config["base_dir"] + "/prokka/{prokka_accession}"
    threads:
        6
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir} --prefix {params.prefix} --cpus {threads} {input} --force
        """
#
#rule prepare_taxonomy_files:
#    input:
#        config["base_dir"] + "/" + config["prefix"] + ".sample_gtdb.metadata.tsv"
#    output:
#        names=config["base_dir"] + "/taxonomy_data/names.dmp",
#        nodes=config["base_dir"] + "/taxonomy_data/nodes.dmp",
#        taxonmap=config["base_dir"] + "/taxonomy_data/taxid.map"
#    shell:
#        """
#        python scripts/create_taxon_data.py --metadata {input} \
#            --out-nodes {output.nodes} \
#            --out-names {output.names} \
#            --out-taxid {output.taxonmap}
#        """
#
#rule create_protein_to_taxa_map:
#    input:
#        taxonmap=config["base_dir"] + "/taxonomy_data/taxid.map",
#        proteomes=config["base_dir"] + "/" + config["prefix"] + ".sample_gtdb.proteomes/"
#    output:
#        config["base_dir"] + "/taxonomy_data/prot2taxid.map"
#    conda:
#        "../envs/biopython.yaml"
#    shell:
#        """
#        python scripts/create_taxonmap.py {input.taxonmap} {input.proteomes} {output}
#        """
#

rule make_diamond_db:
    input:
        prokka=get_prokka_proteomes,
        ncbi=get_ncbi_proteomes
        #names=config["base_dir"] + "/taxonomy_data/names.dmp",
        #nodes=config["base_dir"] + "/taxonomy_data/nodes.dmp",
        #taxonmap=config["base_dir"] + "/taxonomy_data/prot2taxid.map"
    output:
        config["base_dir"] + "/diamond_db/{prefix}.{method}.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    shell:
        """
        cat {input.prokka} {input.ncbi} | diamond makedb --db {output} -p {threads}
        """
        #--taxonnames {input.names} --taxonnodes {input.nodes} --taxonmap {input.taxonmap}


rule make_blast_db:
    input:
        config["base_dir"] + "/{prefix}.{method}.proteomes/"
    output:
        config["base_dir"] + "/ncbi_blastp_db/{prefix}.{method}.pdb"
    params:
        prefix=config["base_dir"] + "/ncbi_blastp_db/{prefix}.{method}",
        title="{prefix}.{method}"
    conda:
        "../envs/ncbi_blast.yaml"
    shell:
        """
        cat {input}/* | makeblastdb -in - -dbtype prot -out {params.prefix} -title {params.title};
        """
