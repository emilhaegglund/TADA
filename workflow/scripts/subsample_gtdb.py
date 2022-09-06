import pandas as pd
import argparse


def read_command_line():
    """Read command line arguments"""
    args = argparse.ArgumentParser()
    args.add_argument("--gtdb-taxonomy", required=True)
    args.add_argument("--gtdb-metadata", required=True)
    args.add_argument("--output-refseq", required=True)
    args.add_argument("--output-genbank", required=True)
    args.add_argument(
        "--taxonomic-level",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        default="class",
    )
    args.add_argument("--max-taxa", type=int, default=1)
    args.add_argument("--completeness", type=float, default=0)
    args.add_argument("--contamination", type=float, default=100)
    args.add_argument("--gtdb-representative")

    return args.parse_args()


args = read_command_line()

# Read GTDB taxonomy
taxa_df = pd.read_csv(
    args.gtdb_taxonomy, sep="\t", names=["assembly_accession", "taxonomy"]
)
taxa_df[
    ["domain", "phylum", "class", "order", "family", "genus", "species"]
] = taxa_df.taxonomy.str.split(";", expand=True)

# Read GTDB metadata
metadata_df = pd.read_csv(args.gtdb_metadata, sep="\t")

# Merge the tables
df = pd.merge(
    left=taxa_df,
    right=metadata_df,
    left_on="assembly_accession",
    right_on="accession",
)

print(df.shape)
# Filter based on contamination and completeness
df = df[
    (df.checkm_contamination <= args.contamination)
    & (df.checkm_completeness >= args.completeness)
]

if args.gtdb_representative == True:
    df = df[df.gtdb_representative == "t"]

# Remove accession prefix
df["accession"] = df["accession"].str.replace("GB_", "")
df["accession"] = df["accession"].str.replace("RS_", "")


sampled_accessions = []
for i, taxa_level_df in df.groupby(args.taxonomic_level):
    # Can't take a sample if the sample size we ask for is larger than
    # the number of taxa in that group. In that case, use all taxa in
    # the group.
    if taxa_level_df.shape[0] > args.max_taxa:
        sampled_accessions.append(taxa_level_df.sample(args.max_taxa))
    else:
        sampled_accessions.append(taxa_level_df)


sampled_df = pd.concat(sampled_accessions)
sampled_df_genbank = sampled_df[sampled_df["accession"].str.contains("GCA")]
sampled_df_refseq = sampled_df[sampled_df["accession"].str.contains("GCF")]
sampled_df_refseq["accession"].to_csv(args.output_refseq, sep="\t", index=False)
sampled_df_genbank["accession"].to_csv(args.output_genbank, sep="\t", index=False)
