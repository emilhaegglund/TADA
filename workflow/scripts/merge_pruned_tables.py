import pandas as pd

bac120_df = pd.read_csv(snakemake.input.bacteria_metadata, sep="\t")
ar53_df = pd.read_csv(snakemake.input.archaea_metadata, sep="\t")

df = pd.concat([bac120_df, ar53_df])
df.to_csv(snakemake.output.metadata, sep="\t", index=False)