rule download_gtdb_metadata:
    output:
        "gtdb_data/{domain}_metadata_r{version}.raw.tsv"
    params:
        output_dir="gtdb_data/",
        base_url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release{version}/{version}.0/"
        #version=config["sample_gtdb"]["version"]
    conda:
        "../envs/base.yaml"
    log:
        "logs/dowload_gtdb_metadata.{domain}_{version}.log"
    shell:
        """
        if [[ {wildcards.version} == "220" ]]; then
            filename="{wildcards.domain}_metadata_r{wildcards.version}.tsv.gz"
        else
            filename="{wildcards.domain}_metadata_r{wildcards.version}.tar.gz"
        fi
        
        url="{params.base_url}$filename"

        wget -P {params.output_dir} $url && \
            gunzip -c {params.output_dir}/$filename > {output} 2> {log}
        """

rule download_gtdb_phylogenies:
    output:
       "gtdb_data/{domain}_r{version}.tree"
    params:
        output_dir="gtdb_data/",
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release{version}/{version}.0/{domain}_r{version}.tree"
    conda:
        "../envs/base.yaml"
    log:
        "logs/dowload_gtdb_phylogenies.{domain}_{version}.log"
    shell:
        """
        wget -P {params.output_dir} {params.url} 2> {log}
        """

rule remove_suppressed_records_from_metadata:
    input:
        metadata="gtdb_data/{domain}_metadata_r{version}.raw.tsv",
        suppressed_genbank_records="ncbi_data/assembly_summary_genbank_historical.txt",
        suppressed_refseq_records="ncbi_data/assembly_summary_refseq_historical.txt"
    output:
        "gtdb_data/{domain}_metadata_r{version}.wo_suppressed_records.tsv",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/remove_suppressed_records_from_metadata.py"

rule merge_metadata_tables:
    input:
        bac_metadata="gtdb_data/bac120_metadata_r{version}.wo_suppressed_records.tsv",
        ar_metadata="gtdb_data/ar53_metadata_r{version}.wo_suppressed_records.tsv"
    output:
        "gtdb_data/metadata_r{version}.wo_suppressed_records.tsv"
    conda:
        "../envs/biopython.yaml"
    log:
        "logs/merge_metadata_tables.{version}.log"
    script:
        "../scripts/merge_gtdb_ar_bac_metadata.py" 