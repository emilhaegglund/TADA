rule download_assembly_summary:
    output:
        config["paths"]["results"] + "/refseq_data/assembly_summary_refseq.txt"
    params:
        outdir=config["paths"]["results"] + "/refseq_data/assembly_summary_refseq.txt"
    shell:
        "wget -P {params.outdir} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
