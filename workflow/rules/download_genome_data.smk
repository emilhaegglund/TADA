rule get_sampled_accessions:
    input:
       config["base_dir"] + "/{prefix}.{method}.metadata.tsv"
    output:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt",
    shell:
        "cut -f1 {input} | grep -v 'Accessions' > {output}"

rule get_refseq_accessions:
    input:
        config["base_dir"] + "/{prefix}.{method}.sampled_accessions.txt",
    output:
        temporary(config["base_dir"] + "/{prefix}.{method}.refseq.accessions.txt")
    shell:
        """
        echo 'accession' > {output};
        grep 'GCF' {input} >> {output} || true;
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

rule ncbi_download_proteomes:
    input:
        config["base_dir"] + "/{prefix}.{method}.{database}.accessions.txt"
    output:
        directory(config["base_dir"] + "/{prefix}.{method}.{database}.proteome_data")
    params:
        db='{database}',
        output_dir=config["base_dir"] + "/{prefix}.{method}.{database}.proteome_data"
    threads:
        12
    conda:
        "../envs/ncbi-download-genomes.yaml"
    shell:
        """
        mkdir -p {output}
        ncbi-genome-download -F protein-fasta -A {input} --flat-output -p {threads} -o {output} -s {params.db} all || true
        """

# Not all assemblies in GenBank contains annotation. Find the sampled taxa for which it was not possible
# to download protein-sequences so we can annotate them with prodigal.
rule  get_accessions_wo_annotation:
    input:
        proteome_dir = config["base_dir"] + "/" + config["prefix"] + "." + config["method"] + ".genbank.proteome_data",
        genbank_accessions = config["base_dir"] + "/" + config["prefix"] + "." + config["method"] + ".genbank.accessions.txt"
    output:
        config["base_dir"] + "/genbank_genomes_wo_annotation.tsv"
    shell:
        """
        python scripts/check_proteome_download.py {input.genbank_accessions} {input.proteome_dir} {output}
        """

# Needs to define a checkpoint here since we don't know for which taxa we have to run prodigal for.
checkpoint ncbi_download_genomes:
    input:
        config["base_dir"] + "/genbank_genomes_wo_annotation.tsv"
    output:
        directory(config["base_dir"] + "/genbank_genomes/")
    params:
        db='genbank'
    threads:
        12
    conda:
        "../envs/ncbi-download-genomes.yaml"
    shell:
        """
        mkdir -p {output}
        ncbi-genome-download -F fasta -A {input} --flat-output -p {threads} -o {output} -s {params.db} all || true
        """

# Input function
def get_genome_accessions(wildcards):
    ck_output = checkpoints.ncbi_download_genomes.get(**wildcards).output[0]
    return expand(config["base_dir"] + "/{annotation}/{accession}/{accession}.faa.gz",
                    annotation='Prokka',
                    accession=glob_wildcards(os.path.join(ck_output, "{accession}_genomic.fna.gz")).accession)

def get_accession(wildcards):
    return "_".join(wildcards.accession.split('_'))[:2]

rule unzip_genomes:
    input:
        config["base_dir"] + "/genbank_genomes/{accession}_genomic.fna.gz"
    output:
        temp(config["base_dir"] + "/genbank_genomes_unzipped/{accession}_genomic.fna")
    shell:
        "gunzip -c {input} > {output}"

rule prokka:
    input:
        fasta=config["base_dir"] + "/genbank_genomes_unzipped/{accession}_genomic.fna",
        metadata=config["base_dir"] + "/" + config["prefix"] + "." + config["method"] + ".metadata.tsv"
    output:
        config["base_dir"] + "/Prokka/{accession}/{accession}.faa"
    params:
        short_accession=lambda wc: '_'.join(wc.get("accession").split('_')[0:2]),
        accession="{accession}",
        outdir=config["base_dir"] + "/Prokka/{accession}"
    threads:
        6
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        echo {params.short_accession};
        KINGDOM=$(awk -F'\\t' '$1=="{params.short_accession}" {{ print $112 }}' {input.metadata} | uniq);
        echo $KINGDOM;
        prokka --kingdom $KINGDOM --outdir {params.outdir} --prefix {params.accession} --cpus {threads} {input} --force
        """

rule zip_prokka:
    input:
        config["base_dir"] + "/Prokka/{accession}/{accession}.faa"
    output:
        config["base_dir"] + "/Prokka/{accession}/{accession}.faa.gz"
    shell:
        "gzip {input}"

rule gather_protein_sequences:
    input:
        annotation=get_genome_accessions,
        refseq=config["base_dir"] + "/{prefix}.{method}.refseq.proteome_data/",
        genbank=config["base_dir"] + "/{prefix}.{method}.genbank.proteome_data/"
    output:
        directory(config["base_dir"] + "/{prefix}.{method}.proteomes/")
    shell:
        """
        mkdir {output};
        mv {input.refseq}/* {input.genbank}/* {output} || true;
        cp {input.annotation} {output} || true;
        """

rule prepare_taxonomy_files:
    input:
        config["base_dir"] + "/" + config["prefix"] + ".sample_gtdb.metadata.tsv"
    output:
        names=config["base_dir"] + "/taxonomy_data/names.dmp",
        nodes=config["base_dir"] + "/taxonomy_data/nodes.dmp",
        taxonmap=config["base_dir"] + "/taxonomy_data/taxid.map"
    shell:
        """
        python scripts/create_taxon_data.py --metadata {input} \
            --out-nodes {output.nodes} \
            --out-names {output.names} \
            --out-taxid {output.taxonmap}
        """

rule create_protein_to_taxa_map:
    input:
        taxonmap=config["base_dir"] + "/taxonomy_data/taxid.map",
        proteomes=config["base_dir"] + "/" + config["prefix"] + ".sample_gtdb.proteomes/"
    output:
        config["base_dir"] + "/taxonomy_data/prot2taxid.map"
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python scripts/create_taxonmap.py {input.taxonmap} {input.proteomes} {output}
        """

rule make_diamond_db:
    input:
        proteomes=config["base_dir"] + "/{prefix}.{method}.proteomes/",
        #names=config["base_dir"] + "/taxonomy_data/names.dmp",
        #nodes=config["base_dir"] + "/taxonomy_data/nodes.dmp",
        #taxonmap=config["base_dir"] + "/taxonomy_data/prot2taxid.map"
    output:
        config["base_dir"] + "/{prefix}.{method}.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    shell:
        """
        zcat {input.proteomes}/* | diamond makedb --db {output} -p {threads}
        """
        #--taxonnames {input.names} --taxonnodes {input.nodes} --taxonmap {input.taxonmap}


rule make_blast_db:
    input:
        config["base_dir"] + "/{prefix}_{method}_proteomes/"
    output:
        config["base_dir"] + "/{prefix}.{method}.pdb"
    params:
        prefix=config["base_dir"] + "/{prefix}.{method}",
        title="{prefix}.{method}"
    conda:
        "../envs/ncbi_blast.yaml"
    shell:
        """
        zcat {input}/* | makeblastdb -in - -dbtype prot -out {params.prefix} -title {params.title};
        """
#rule make_mmseqs_db:
#    input:
#
#    output
#
#    params:
#
#    conda:
#
#    threads:
#
#    shell:
