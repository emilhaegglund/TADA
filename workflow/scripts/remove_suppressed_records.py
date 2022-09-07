import pandas as pd
import argparse


def parse_command_line():
    args = argparse.ArgumentParser()
    args.add_argument("--metadata", required=True)
    args.add_argument("--suppressed-records", required=True)
    args.add_argument("--output", required=True)

    return args.parse_args()


args = parse_command_line()
metadata_df = pd.read_csv(args.metadata, sep="\t", low_memory=False)

suppressed_names = [
    "assembly_accession",
]

suppressed_df = pd.read_csv(
    args.suppressed_records, sep="\t", comment="#", names=suppressed_names, usecols=[0]
)

suppressed_accessions = suppressed_df["assembly_accession"].to_list()

metadata_df = metadata_df[
    ~metadata_df["ncbi_genbank_assembly_accession"].isin(suppressed_accessions)
]
metadata_df.to_csv(args.output, sep="\t", index=False)
