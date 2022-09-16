rule download_gtdb_metadata:
    output:
       config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.tsv"
    params:
        output_dir=config["paths"]["results"] + "/gtdb_data/",
        zip_file=config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.tar.gz",
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/{domain}_metadata_r207.tar.gz"
    shell:
        """
        wget -P {params.output_dir} {params.url} && \
        tar -xvf {params.zip_file} -C {params.output_dir}
        """

rule download_gtdb_phylogenies:
    output:
       config["paths"]["results"] + "/gtdb_data/{domain}_r207.tree"
    params:
        output_dir=config["paths"]["results"] + "/gtdb_data/",
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/{domain}_r207.tree"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule download_gtdb_taxdump:
    output:
        directory(config["paths"]["results"] + "/gtdb_data/gtdb-taxdump/")
    params:
        output_dir=config["paths"]["results"] + "/gtdb_data/",
        zip_file=config["paths"]["results"] + "/gtdb_data/gtdb-taxdump.tar.gz",
        url="https://github.com/shenwei356/gtdb-taxdump/releases/download/v0.1.1/gtdb-taxdump.tar.gz"
    shell:
        """
        wget -P {params.output_dir} {params.url} && \
        tar -xvf {params.zip_file} -C {params.output_dir}
        """

rule download_suppressed_genbank_records:
    output:
        config["paths"]["results"] + "/gtdb_data/assembly_summary_genbank_historical.txt"
    params:
        output_dir=config["paths"]["results"] + "/gtdb_data/",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank_historical.txt"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule download_suppressed_refseq_records:
    output:
        config["paths"]["results"] + "/gtdb_data/assembly_summary_refseq_historical.txt"
    params:
        output_dir=config["paths"]["results"] + "/gtdb_data/",
        url="https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq_historical.txt"
    shell:
        """
        wget -P {params.output_dir} {params.url}
        """

rule remove_suppressed_records_from_metadata:
    input:
        metadata=config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.tsv",
        suppressed_genbank_records=config["paths"]["results"] + "/gtdb_data/assembly_summary_genbank_historical.txt",
        suppressed_refseq_records=config["paths"]["results"] + "/gtdb_data/assembly_summary_refseq_historical.txt"
    output:
        config["paths"]["results"] + "/gtdb_data/{domain}_metadata_r207.wo_suppressed_records.tsv",
    shell:
        """
        python scripts/remove_suppressed_records_from_metadata.py \
            --metadata {input.metadata} \
            --suppressed-genbank-records {input.suppressed_genbank_records} \
            --suppressed-refseq-records {input.suppressed_refseq_records} \
            --output {output}
        """

rule merge_metadata_tables:
    input:
        bac_metadata=config["paths"]["results"] + "/gtdb_data/bac120_metadata_r207.wo_suppressed_records.tsv",
        ar_metadata=config["paths"]["results"] + "/gtdb_data/ar53_metadata_r207.wo_suppressed_records.tsv"
    output:
        config["paths"]["results"] + "/gtdb_data/metadata_r207.wo_suppressed_records.tsv"
    shell:
        """
        python scripts/merge_tables.py {input.bac_metadata} {input.ar_metadata} {output}
        """
