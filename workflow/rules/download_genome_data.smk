rule get_refseq_accessions:
    input:
        config["paths"]["results"] + "/{method}/sampled_accessions.txt",
    output:
        config["paths"]["results"] + "/{method}/sampled_refseq_accessions.txt"
    shell:
        """
        echo 'accession' > {output};
        grep 'GCF' {input} >> {output} || true;
        """

rule get_genbank_accessions:
    input:
        config["paths"]["results"] + "/{method}/sampled_accessions.txt",
    output:
        temporary(config["paths"]["results"] + "/{method}/sampled_genbank_accessions.txt")
    shell:
        """
        echo 'accession' > {output};
        grep 'GCA' {input} >> {output} || true;
        """

rule ncbi_download_proteomes:
    input:
        config["paths"]["results"] + "/{method}/sampled_{database}_accessions.txt"
    output:
        directory(config["paths"]["results"] + "/{method}/sampled_{database}_proteomes/")
    params:
        db='{database}',
        output_dir=config["paths"]["results"] + "/{method}/sampled_{database}_proteomes/"
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
        proteome_dir = rules.ncbi_download_proteomes.output,
        genbank_accessions = config["paths"]["results"] + "/{method}/sampled_{database}_accessions.txt"
    output:
        config["paths"]["results"] + "/{method}/sampled_{database}_genomes_wo_annotation.tsv"
    shell:
        """
        python scripts/check_proteome_download.py {input.genbank_accessions} {input.proteome_dir} {output}
        """

# Needs to define a checkpoint here since we don't know for which taxa we have to run prodigal for.
checkpoint ncbi_download_genomes:
    input:
        config["paths"]["results"] + "/{method}/sampled_genbank_genomes_wo_annotation.tsv"
    output:
        directory(config["paths"]["results"] + "/{method}/sampled_genbank_genomes/")
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
    return expand(config["paths"]["results"] + "/{method}/{annotation}/{accession}/{accession}.faa",
                    method=wildcards.method,
                    annotation=config["annotation_method"],
                    accession=glob_wildcards(os.path.join(ck_output, "{accession}_genomic.fna.gz")).accession)

def get_accession(wildcards):
    return "_".join(wildcards.accession.split('_'))[:2]

rule unzip_genomes:
    input:
        config["paths"]["results"] + "/{method}/sampled_genbank_genomes/{accession}_genomic.fna.gz"
    output:
        temp(config["paths"]["results"] + "/{method}/sampled_genbank_genomes_unzipped/{accession}_genomic.fna")
    shell:
        "gunzip -c {input} > {output}"

#rule prodigal:
#    input:
#        config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes_unzipped/{accession}_genomic.fna"
#    output:
#        config["paths"]["results"] + "/{method}/prodigal_{domain}/{accession}.faa"
#    threads:
#        1
#    conda:
#        "../envs/prodigal.yaml"
#    shell:
#        """
#        prodigal -i {input} -a {output} -q
#        """

rule prokka:
    input:
        fasta=config["paths"]["results"] + "/{method}/sampled_genbank_genomes_unzipped/{accession}_genomic.fna",
        metadata=config["paths"]["results"] + "/{method}/sampled_accessions.metadata.tsv"
    output:
        config["paths"]["results"] + "/{method}/Prokka/{accession}/{accession}.faa"
    params:
        accession=lambda wc: '_'.join(wc.get("accession").split('_')[0:2]),
        full_accession="{accession}",
        outdir=config["paths"]["results"] + "/{method}/Prokka/{accession}"
    threads:
        6
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        echo {params.accession};
        KINGDOM=$(awk -F'\\t' '$10=="{params.accession}" {{ print $3 }}' {input.metadata} | uniq);
        echo $KINGDOM;
        prokka --kingdom $KINGDOM --outdir {params.outdir} --prefix {params.full_accession} --cpus {threads} {input} --force
        """

#rule bakta:
#    input:
#        config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes_unzipped/{accession}_genomic.fna"
#    output:
#        config["paths"]["results"] + "/{method}/bakta_{domain}/{accession}.faa"
#    threads:
#        6
#    conda:
#        "../envs/bakta.yaml"


rule gather_protein_sequences:
    input:
        refseq=config["paths"]["results"] + "/{method}/sampled_refseq_proteomes/",
        genbank=config["paths"]["results"] + "/{method}/sampled_genbank_proteomes/"
    output:
        temp(config["paths"]["results"] + "/{method}/sampled_proteomes.faa")
    shell:
        """
        zcat {input.refseq}/* {input.genbank}/* > {output} || true
        """

rule make_diamond_db:
    input:
        annotated_proteomes=get_genome_accessions,
        downloaded_proteomes=config["paths"]["results"] + "/{method}/sampled_proteomes.faa",
    output:
        config["paths"]["results"] + "/{method}/{method}.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    shell:
        """
        cat  {input.annotated_proteomes} {input.downloaded_proteomes} | diamond makedb --db {output} -p {threads}
        """
