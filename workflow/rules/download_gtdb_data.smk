rule download_gtdb_taxonomy:
    output:
        config["paths"]["results"] + "/gtdb_data/{domain}_taxonomy_r207.tsv"
    params:
        output_dir=config["paths"]["results"] + "/gtdb_data/",
        zip_file=config["paths"]["results"] + "/gtdb_data/{domain}_taxonomy_r207.tsv.gz",
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/{domain}_taxonomy_r207.tsv.gz"
    shell:
        """
        wget -P {params.output_dir} {params.url} && \
        gunzip {params.zip_file}
        """

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
