import pandas as pd
import sys

dfs = []
for df_path in sys.argv[1:-1]:
    print(df_path)
    df = pd.read_csv(df_path, sep='\t')
    dfs.append(df)

df = pd.concat(dfs)
df = df[["Assembly Accession", "Organism Taxonomic ID"]]

df.drop_duplicates(inplace=True)
df.rename(columns={"Assembly Accession":"assembly_accession", "Organism Taxonomic ID":"taxid"}, inplace=True)
print(df.columns)
print(df)
print(sys.argv[-1])
# Select only RefSeq
df = df[df["assembly_accession"].str.contains("GCF")]
df.to_csv(sys.argv[-1], sep='\t', index=False)
