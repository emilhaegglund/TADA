
rule subsample_assembly_summary:
    input:
        metadata="ncbi_data/datasets.tsv",
        names="ncbi_data/taxdmp/names.dmp",
        nodes="ncbi_data/taxdmp/nodes.dmp",
    output:
        "sample_ncbi.annotation_data.tsv"
    params:
        sampling_scheme=config["sample_ncbi"]["sampling_scheme"]
    script:
        "../scripts/subsample_refseq.py"
