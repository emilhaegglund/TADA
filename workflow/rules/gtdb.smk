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
        taxonomy=config["paths"]["results"] + "/gtdb_data/{domain}_taxonomy_r207.wo_suppressed_records.tsv",
        metadata=config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.wo_suppressed_records.tsv"
    output:
        genbank=config["paths"]["results"] + "/subsample_gtdb/sampled_{domain}_genbank_accessions.tsv",
        refseq=config["paths"]["results"] + "/subsample_gtdb/sampled_{domain}_refseq_accessions.tsv"
    params:
        taxonomic_level=config["subsample_gtdb"]["taxonomic_level"],
        max_taxa=config["subsample_gtdb"]["max_taxa"],
        completeness=config["subsample_gtdb"]["completeness"],
        contamination=config["subsample_gtdb"]["contamination"],
        gtdb_representative=config["subsample_gtdb"]["gtdb_representative"]
    shell:
        """
        python scripts/subsample_gtdb.py \
            --gtdb-taxonomy {input.taxonomy} \
            --gtdb-metadata {input.metadata} \
            --taxonomic-level {params.taxonomic_level} \
            --max-taxa {params.max_taxa} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            --gtdb-representative {params.gtdb_representative} \
            --output-refseq {output.refseq} \
            --output-genbank {output.genbank}
        """

rule ncbi_download_proteomes:
    input:
        config["paths"]["results"] + "/{method}/sampled_{domain}_{database}_accessions.tsv"
    output:
        directory(config["paths"]["results"] + "/{method}/sampled_{domain}_{database}_proteomes/")
    params:
        db='{database}',
        output_dir=config["paths"]["results"] + "/{method}/sampled_{domain}_{database}_proteomes/"
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
        genbank_accessions = config["paths"]["results"] + "/{method}/sampled_{domain}_{database}_accessions.tsv"
    output:
        config["paths"]["results"] + "/{method}/sampled_{domain}_{database}_genomes_wo_annotation.tsv"
    shell:
        """
        python scripts/check_proteome_download.py {input.genbank_accessions} {input.proteome_dir} {output}
        """

# Needs to define a checkpoint here since we don't know for which taxa we have to run prodigal for.
checkpoint ncbi_download_genomes:
    input:
        config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes_wo_annotation.tsv"
    output:
        directory(config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes/")
    params:
        db='genbank'
    threads:
        12
    conda:
        "../envs/ncbi-download-genomes.yaml"
    shell:
        """
        ncbi-genome-download -F fasta -A {input} --flat-output -p {threads} -o {output} -s {params.db} all
        """

# Input function
def get_genome_accessions(wildcards):
    ck_output = checkpoints.ncbi_download_genomes.get(**wildcards).output[0]
    return expand(config["paths"]["results"] + "/{method}/prodigal_{domain}/{accession}.faa",
                    method=wildcards.method,
                    domain=wildcards.domain,
                    accession=glob_wildcards(os.path.join(ck_output, "{accession}_genomic.fna.gz")).accession)

rule unzip_genomes:
    input:
        config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes/{accession}_genomic.fna.gz"
    output:
        temp(config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes_unzipped/{accession}_genomic.fna")
    shell:
        "gunzip -c {input} > {output}"

rule prodigal:
    input:
        config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_genomes_unzipped/{accession}_genomic.fna"
    output:
        config["paths"]["results"] + "/{method}/prodigal_{domain}/{accession}.faa"
    threads:
        1
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input} -a {output} -q
        """

rule gather_prodigal_protein_sequences:
    input:
        get_genome_accessions
    output:
        temp(config["paths"]["results"] + '/{method}/{domain}_prodigal.faa')
    shell:
        "cat {input} > {output}"

rule gather_protein_sequences:
    input:
        refseq=config["paths"]["results"] + "/{method}/sampled_{domain}_refseq_proteomes/",
        genbank=config["paths"]["results"] + "/{method}/sampled_{domain}_genbank_proteomes/"
    output:
        temp(config["paths"]["results"] + "/{method}/sampled_{domain}.faa")
    shell:
        """
        zcat {input.refseq}/* {input.genbank}/* > {output} || true
        """

rule make_diamond_db:
    input:
        archaea_prodigal=config["paths"]["results"] + "/{method}/ar53_prodigal.faa",
        bacteria_prodigal=config["paths"]["results"] + "/{method}/bac120_prodigal.faa",
        bacteria=config["paths"]["results"] + "/{method}/sampled_bac120.faa",
        archaea=config["paths"]["results"] + "/{method}/sampled_ar53.faa"
    output:
        config["paths"]["results"] + "/{method}/{method}.dmnd"
    conda:
        "../envs/diamond.yaml"
    threads:
        12
    shell:
        """
        cat {input.archaea} {input.bacteria} {input.archaea_prodigal} {input.bacteria_prodigal} | diamond makedb --db {output} -p {threads}
        """
