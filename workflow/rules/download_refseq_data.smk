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
