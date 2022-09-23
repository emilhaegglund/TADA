
rule subsample_assembly_summary:
    input:
        metadata=config["base_dir"] + "/ncbi_data/assembly_summary_refseq.txt",
        assembly_summary_historical=config["base_dir"] + "/ncbi_data/assembly_summary_refseq_historical.txt",
        names=config["base_dir"] + "/ncbi_data/taxdmp/names.dmp",
        nodes=config["base_dir"] + "/ncbi_data/taxdmp/nodes.dmp",
    output:
        config["base_dir"] + "/{prefix}.sample_refseq.metadata.tsv"
    params:
        sampling_scheme=config["sample_refseq"]["sampling_scheme"]
    shell:
        """
        python scripts/subsample_refseq.py \
            --refseq-metadata {input.metadata} \
            --sampling-scheme {params.sampling_scheme} \
            --names {input.names} \
            --nodes {input.nodes} \
            --historical {input.assembly_summary_historical} \
            --output {output}
        """
