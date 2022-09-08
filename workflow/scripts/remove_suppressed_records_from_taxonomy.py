import pandas as pd
import argparse


def parse_command_line():
    args = argparse.ArgumentParser()
    args.add_argument("--taxonomy", required=True)
    args.add_argument("--suppressed-genbank-records", required=True)
    args.add_argument("--suppressed-refseq-records", required=True)
    args.add_argument("--output", required=True)

    return args.parse_args()


args = parse_command_line()
taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", low_memory=False, names=["accession", "taxonomy"])

suppressed_names = [
    "assembly_accession",
]

suppressed_refseq_df = pd.read_csv(
    args.suppressed_refseq_records,
    sep="\t",
    comment="#",
    names=suppressed_names,
    usecols=[0],
)

suppressed_genbank_df = pd.read_csv(
    args.suppressed_genbank_records,
    sep="\t",
    comment="#",
    names=suppressed_names,
    usecols=[0],
)

# For suppressed RefSeq records, fall back to genbank-accessions
taxonomy_df["assembly_accession"] = taxonomy_df["accession"]
taxonomy_df["assembly_accession"] = taxonomy_df["assembly_accession"].str.replace(
    "RS_", ""
)
taxonomy_df["assembly_accession"] = taxonomy_df["assembly_accession"].str.replace(
    "GB_", ""
)
suppressed_refseq_accessions = suppressed_refseq_df["assembly_accession"].to_list()
taxonomy_df.loc[
    (taxonomy_df["assembly_accession"].isin(suppressed_refseq_accessions)), "accession"
] = taxonomy_df["accession"].str.replace("RS_GCF", "GB_GCA")

suppressed_genbank_accessions = suppressed_genbank_df["assembly_accession"].to_list()

taxonomy_df = taxonomy_df[
    ~taxonomy_df["assembly_accession"].isin(suppressed_genbank_accessions)
]
taxonomy_df[["accession", "taxonomy"]].to_csv(args.output, sep="\t", index=False)

