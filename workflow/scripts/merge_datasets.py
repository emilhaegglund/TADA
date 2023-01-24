import pandas as pd
import numpy as np
import sys

dfs = []
for df_path in snakemake.input:
    print(df_path)
    df = pd.read_csv(df_path, sep="\t")
    dfs.append(df)

df = pd.concat(dfs)
df = df[df["Assembly Status"] != "suppressed"]
df = df[["Assembly Accession", "Organism Taxonomic ID", "Annotation Name"]]
df["annotation"] = pd.notnull(df['Annotation Name'])
#df = df[df["miss_annotation"] == True]
df.rename(
    columns={
        "Assembly Accession": "assembly_accession",
        "Organism Taxonomic ID": "taxid",
    },
    inplace=True,
)
df.drop_duplicates(subset="assembly_accession", inplace=True)
df[["assembly_accession", "taxid", "annotation"]].to_csv(
    snakemake.output[0], sep="\t", index=False
)
