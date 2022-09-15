import pandas as pd
import sys

bac120_df = pd.read_csv(sys.argv[1], sep='\t')
ar53_df = pd.read_csv(sys.argv[2], sep='\t')

df = pd.concat([bac120_df, ar53_df])
df.to_csv(sys.argv[3], sep='\t', index=False)
