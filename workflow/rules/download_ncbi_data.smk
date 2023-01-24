TAXA = []
if config["method"] == "sample_ncbi":
    with open(config["sample_ncbi"]["sampling_scheme"], "r") as stream:
        sampling_scheme = yaml.safe_load(stream)
    for taxa in sampling_scheme.keys():
        TAXA.append(taxa.replace(" ", "_"))

rule download_summary:
    """
    Use t
    """
    output:
        "ncbi_data/{taxa}.tsv"
    params:
        source=config["sample_ncbi"]["source"]
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        taxa=$(echo {wildcards.taxa} | sed -e "s/_/ /");
        taxa_new="'$taxa'";
        datasets summary genome taxon "$taxa_new" --assembly-source {params.source} \
            --as-json-lines | \
        dataformat tsv genome --fields accession,annotinfo-name,organism-tax-id > {output};
        """

rule merge_genome_summary:
    """
    Merge the genome summary files into a single tsv table
    """
    input:
        expand("ncbi_data/{taxa}.tsv", taxa=TAXA)
    output:
        "ncbi_data/datasets.tsv"
    script:
        "../scripts/merge_datasets.py"


rule download_suppressed_genbank_records:
    output:
        "ncbi_data/assembly_summary_genbank_historical.txt"
    params:
        output_dir="ncbi_data/",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule download_suppressed_refseq_records:
    output:
        "ncbi_data/assembly_summary_refseq_historical.txt"
    params:
        output_dir="ncbi_data/",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule download_taxdmp:
    output:
        "ncbi_data/taxdmp/nodes.dmp",
        "ncbi_data/taxdmp/names.dmp"
    params:
        output_dir="ncbi_data/"
    shell:
        "wget -P {params.output_dir} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip && unzip {params.output_dir}/taxdmp.zip -d {params.output_dir}/taxdmp"
