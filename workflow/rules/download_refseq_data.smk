TAXA = []
if config["method"] == "sample_refseq":
    with open(config["sample_refseq"]["sampling_scheme"], "r") as stream:
        sampling_scheme = yaml.safe_load(stream)
    for taxa in sampling_scheme.keys():
        TAXA.append(taxa.replace(" ", "_"))
rule download_summary:
    """
    Use t
    """
    output:
        config["base_dir"] + "/ncbi_data/{taxa}.tsv"
    conda:
        "../envs/ncbi-datasets.yaml"
    shell:
        """
        taxa=$(echo {wildcards.taxa} | sed -e "s/_/ /");
        taxa_new="'$taxa'";
        datasets summary genome taxon "$taxa_new" --as-json-lines | \
        dataformat tsv genome > {output};
        """

rule merge_genome_summary:
    """
    Merge the genome summary files into a single tsv table
    """
    input:
        expand(config["base_dir"] + "/ncbi_data/{taxa}.tsv", taxa=TAXA)
    output:
        config["base_dir"] + "/ncbi_data/datasets.tsv"
    shell:
        "python scripts/merge_datasets.py {input} {output}"


rule download_assembly_summary:
    output:
        config["base_dir"] + "/ncbi_data/assembly_summary_refseq.txt"
    params:
        outdir=config["base_dir"] + "/ncbi_data/"
    shell:
        "wget -P {params.outdir} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

rule download_suppressed_genbank_records:
    output:
        config["base_dir"] + "/ncbi_data/assembly_summary_genbank_historical.txt"
    params:
        output_dir=config["base_dir"] + "/ncbi_data/",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule download_suppressed_refseq_records:
    output:
        config["base_dir"] + "/ncbi_data/assembly_summary_refseq_historical.txt"
    params:
        output_dir=config["base_dir"] + "/ncbi_data/",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule download_taxdmp:
    output:
        config["base_dir"] + "/ncbi_data/taxdmp/nodes.dmp",
        config["base_dir"] + "/ncbi_data/taxdmp/names.dmp"
    params:
        output_dir=config["base_dir"] + "/ncbi_data/"
    shell:
        "wget -P {params.output_dir} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip && unzip {params.output_dir}/taxdmp.zip -d {params.output_dir}/taxdmp"
