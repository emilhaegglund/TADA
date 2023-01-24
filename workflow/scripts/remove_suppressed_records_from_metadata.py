import pandas as pd


metadata = snakemake.input.metadata
suppressed_genbank_records = snakemake.input.suppressed_genbank_records
suppressed_refseq_records = snakemake.input.suppressed_refseq_records
output = snakemake.output[0]

metadata_df = pd.read_csv(metadata, sep="\t", low_memory=False)

suppressed_names = [
    "assembly_accession",
]

suppressed_refseq_df = pd.read_csv(
    suppressed_refseq_records,
    sep="\t",
    comment="#",
    names=suppressed_names,
    usecols=[0],
)

suppressed_genbank_df = pd.read_csv(
    suppressed_genbank_records,
    sep="\t",
    comment="#",
    names=suppressed_names,
    usecols=[0],
)

# For suppressed RefSeq records, fall back to genbank-accessions
metadata_df["assembly_accession"] = metadata_df["accession"]
metadata_df["assembly_accession"] = metadata_df["assembly_accession"].str.replace(
    "RS_", ""
)
metadata_df["assembly_accession"] = metadata_df["assembly_accession"].str.replace(
    "GB_", ""
)
suppressed_refseq_accessions = suppressed_refseq_df["assembly_accession"].to_list()
metadata_df.loc[
    (metadata_df["assembly_accession"].isin(suppressed_refseq_accessions)), "accession"
] = metadata_df["accession"].str.replace("RS_GCF", "GB_GCA")

suppressed_genbank_accessions = suppressed_genbank_df["assembly_accession"].to_list()

metadata_df = metadata_df[
    ~metadata_df["ncbi_genbank_assembly_accession"].isin(suppressed_genbank_accessions)
]
metadata_df.to_csv(output, sep="\t", index=False)
