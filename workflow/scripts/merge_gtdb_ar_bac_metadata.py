import pandas as pd

bac120_df = pd.read_csv(snakemake.input.bac_metadata, sep='\t', low_memory=False)
ar53_df = pd.read_csv(snakemake.input.ar_metadata, sep='\t', low_memory=False)

df = pd.concat([bac120_df, ar53_df])
df.to_csv(snakemake.output[0], sep='\t', index=False)
