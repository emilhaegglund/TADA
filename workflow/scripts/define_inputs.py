import os


def gtdb_data_files(config):
    """
    Define files that needs to be downloaded to run the GTDB-related
    workflows.
    """
    wanted_input = []

    base_dir = config["paths"]["results"]
    gtdb_data_files = [
        "ar53_taxonomy_r207.tsv",
        "bac120_taxonomy_r207.tsv",
        "ar53_metadata_r207.wo_suppressed_records.tsv",
        "bac120_metadata_r207.wo_suppressed_records.tsv",
        "ar53_r207.tree",
        "bac120_r207.tree",
        "assembly_summary_genbank_historical.txt",
        "assembly_summary_refseq_historical.txt",
        "gtdb-taxdump/",
    ]
    for f in gtdb_data_files:
        wanted_input.append(os.path.join(base_dir, "gtdb_data", f))

    return wanted_input


def gtdb_workflow_files(config):
    """
    Define common output files for the subsample and prune workflow.
    """
    wanted_input = []
    base_dir = config["paths"]["results"]
    if config["method"]["subsample_gtdb"]:
        wanted_input.append(
            os.path.join(base_dir, "subsample_gtdb", "subsample_gtdb.dmnd")
        )
    elif config["method"]["prune_gtdb"]:
        wanted_input.append(os.path.join(base_dir, "prune_gtdb", "prune_gtdb.dmnd"))

    return wanted_input
