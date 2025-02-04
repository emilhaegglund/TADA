if "sample_ncbi" not in config.keys():
    config["sample_ncbi"] = {"sampling_scheme": "", "database": "GenBank"}

TAXA = []
if config["method"] == "sample_ncbi":
    with open(config["sample_ncbi"]["sampling_scheme"], "r") as stream:
        sampling_scheme = yaml.safe_load(stream)
    for taxa in sampling_scheme.keys():
        if taxa == 'all':
            TAXA.append("Bacteria")
            TAXA.append("Archaea")
        else:
            TAXA.append(taxa.replace(" ", "_"))

rule download_summary:
    output:
        "ncbi_data/taxa/{taxa}.tsv"
    params:
        source=config["sample_ncbi"]["database"]
    conda:
        "../envs/ncbi-datasets.yaml"
    log:
        "logs/download_summary.{taxa}.log"
    shell:
        """
        taxa=$(echo {wildcards.taxa} | sed -e "s/_/ /");
        datasets summary genome taxon "$taxa" --tax-exact-match --assembly-source {params.source} \
            --as-json-lines | \
        dataformat tsv genome --fields accession,annotinfo-name,assminfo-status,organism-tax-id > {output}; 2> {log}
        """

rule merge_genome_summary:
    """
    Merge the genome summary files into a single tsv table
    """
    input:
        expand("ncbi_data/taxa/{taxa}.tsv", taxa=TAXA)
    output:
        "ncbi_data/datasets_unchecked.tsv"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/merge_datasets.py"

rule check_for_euk:
    input:
        dataset = "ncbi_data/datasets_unchecked.tsv",
        nodes = "ncbi_data/taxdmp/nodes.dmp",
        merged_nodes = "ncbi_data/taxdmp/merged.dmp"
    output:
        dataset = "ncbi_data/datasets.tsv"
    conda:
        "../envs/base.yaml"
    params:
        context = "sampling_scheme"
    script:
        "../scripts/check_for_non_supported_taxa.py"

rule download_suppressed_genbank_records:
    output:
        "ncbi_data/assembly_summary_genbank_historical.txt"
    conda:
        "../envs/base.yaml"
    params:
        output_dir=lambda w, output: os.path.split(output[0])[0],
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt"
    log:
        "logs/download_supressed_genbank_records.log"
    shell:
        """
        wget -P {params.output_dir} {params.url} 2> {log}
        """

rule download_suppressed_refseq_records:
    output:
        "ncbi_data/assembly_summary_refseq_historical.txt"
    params:
        output_dir=lambda w, output: os.path.split(output[0])[0],
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt"
    log:
        "logs/download_suppressed_refseq_records.log"
    shell:
        """
        wget -P {params.output_dir} {params.url} 2> {log}
        """

rule download_taxdmp:
    output:
        "ncbi_data/taxdmp/nodes.dmp",
        "ncbi_data/taxdmp/names.dmp",
        "ncbi_data/taxdmp/merged.dmp"
    conda:
        "../envs/base.yaml"
    params:
        output_dir=lambda w, output: os.path.split(os.path.splitext(output[0])[0])[0],
    log:
        "logs/download_taxdmp.log"
    shell:
        "wget -P {params.output_dir} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip && unzip {params.output_dir}/taxdmp.zip -d {params.output_dir} 2> {log}"

rule get_required_genomes_taxid:
    input:
        config["required"]
    output:
        "ncbi_data/required_genomes_unchecked.tsv"
    conda:
        "../envs/ncbi-datasets.yaml"
    log:
        dataset="logs/get_required_genomes_taxid_dataset.log",
        dataformat="logs/get_required_genomes_taxid_dataformat.log"
    shell:
        """
        datasets summary genome accession --inputfile {input} --as-json-lines 2> {log.dataset} | \
        dataformat tsv genome  --fields accession,annotinfo-name,assminfo-status,organism-tax-id > {output} 2> {log.dataformat}
        """

rule merge_genome_summary_required_genomes:
    """
    Merge the genome summary files into a single tsv table
    """
    input:
        "ncbi_data/required_genomes_unchecked.tsv"
    output:
        "ncbi_data/required_genomes_unchecked_anno_info.tsv"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/merge_datasets.py"

rule check_required_genomes_for_euk:
    input:
        dataset = "ncbi_data/required_genomes_unchecked_anno_info.tsv",
        nodes = "ncbi_data/taxdmp/nodes.dmp",
        merged_nodes = "ncbi_data/taxdmp/merged.dmp"
    output:
        dataset = "ncbi_data/required_genomes_checked.tsv"
    conda:
        "../envs/base.yaml"
    params:
        context="required_genomes"
    script:
        "../scripts/check_for_non_supported_taxa.py"