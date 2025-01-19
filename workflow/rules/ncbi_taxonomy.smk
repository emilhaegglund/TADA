rule subsample_assembly_summary:
    input:
        metadata="ncbi_data/datasets.tsv",
        names="ncbi_data/taxdmp/names.dmp",
        nodes="ncbi_data/taxdmp/nodes.dmp",
        merged="ncbi_data/taxdmp/merged.dmp",
        required_genomes="ncbi_data/required_genomes_checked.tsv" if config["required"] != "" else []
    output:
        "sample_ncbi.annotation_data.tsv"
    conda:
        "../envs/base.yaml"
    params:
        sampling_scheme=config["sample_ncbi"]["sampling_scheme"],
        seed=config["seed"]
    script:
        "../scripts/subsample_ncbi_taxonomy.py"
