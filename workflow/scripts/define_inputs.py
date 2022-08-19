import os


def gtdb_data_files(config):
    wanted_input = []

    base_dir = config["paths"]["results"]
    gtdb_data_files = [
        "ar53_taxonomy_r207.tsv",
        "bac120_taxonomy_r207.tsv",
        "ar53_metadata_r207.tsv",
        "bac120_metadata_r207.tsv",
    ]
    for f in gtdb_data_files:
        wanted_input.append(os.path.join(base_dir, "gtdb_data", f))

    return wanted_input
