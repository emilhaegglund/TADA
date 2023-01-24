import pandas as pd
import os


df = pd.read_csv(snakemake.input[0], sep="\t")

if snakemake.params.status == "without-annotation":
    df = df[df["annotation"] == False]
elif snakemake.params.status == "with-annotation":
    df = df[df["annotation"] == True]

os.mkdir(snakemake.output[0])

for accession in df["assembly_accession"]:
    output_path = os.path.join(snakemake.output[0], accession + ".txt")
    with open(output_path, "w") as fp:
        pass
