"""
Script to merge the genome accession datasets tables downloaded for each taxonomy.
"""
import pandas as pd


# Open and merge the tables
dfs = []
for df_path in snakemake.input:
    df = pd.read_csv(df_path, sep="\t", lineterminator="\n")
    dfs.append(df)
df = pd.concat(dfs)

# Remove the assembly if it is suppressed
df = df[df["Assembly Status"] != "suppressed"]

# Subset tables
df = df[["Assembly Accession", "Organism Taxonomic ID", "Annotation Name"]]

# Create a columns which tells if there is an annotation associated with the assembly
df["annotation"] = pd.notnull(df['Annotation Name'])

# Rename columns
df.rename(
    columns={
        "Assembly Accession": "assembly_accession",
        "Organism Taxonomic ID": "taxid",
    },
    inplace=True,
)

# Drop assemblies that are present multiple times in the table
df.drop_duplicates(subset="assembly_accession", inplace=True)
df[["assembly_accession", "taxid", "annotation"]].to_csv(
    snakemake.output[0], sep="\t", index=False
)
