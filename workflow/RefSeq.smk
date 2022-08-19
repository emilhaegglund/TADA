

rule all:
    input:
        "results/assembly_summary_refseq.txt",
        "results/assembly_summary_refseq_historical.txt",
        "results/taxdmp/nodes.dmp"

""" Download necessary files from NCBI """
rule download_assembly_summary:
    output:
        "results/assembly_summary_refseq.txt"
    shell:
        "wget -P results// https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

rule download_assembly_summary_historical:
    output:
        "results/assembly_summary_refseq_historical.txt"
    shell:
        "wget -P results/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt"

rule download_taxdmp:
    output:
        "results/taxdmp/nodes.dmp"
    shell:
        "wget -P results/ https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip && unzip results/taxdmp.zip -d results/taxdmp"

rule download_prot_to_taxid:
    """ This will download a 7Gb+ file """
    output:
        "results/prot.accession2taxid.gz"
    shell:
        "wget -P results/ https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"

rule subsample_assembly_summary:
    input:
        assembly_summary=
        assembly_summary_historical=
        taxonomy=
    output:
       "results/sampled_accessions.txt"
    params:
        mode="random",
        max_number=1
    shell:
        """"
        python scripts/sample_refseq.py \
            --level genus \
            --assembly_summary {input.assembly_summary} \
            --assembly_summary_historical {input.assembly_summary_historical} \
            --taxonomy {input.taxonomy} \
            --max_number {params.max_number}
            --mode {params.mode}
            --output {output}
        """

rule download_proteomes:
    input:

    output:

    conda:
        "envs/ncbi-download-genomes.yaml"
    shell:
#
#rule build_diamond_database:
#    input:
#
#    output:
#
#    conda:
#        "envs/diamond.yaml"
#    threads:
#
#    shell:
#        "zcat {input} | diamond makedb --db {output} -p {threads}
