import pandas as pd
import sys
import tada

# Read original dataset
df = pd.read_csv(snakemake.input.dataset, sep="\t")
if "Organism Taxonomic ID" in df.columns:
    df.rename(columns={"Organism Taxonomic ID":"taxid"}, inplace=True)

# Read NCBI taxonomy
taxonomy = tada.taxdmp_taxonomy(snakemake.input.nodes)

# Fix merged nodes
merged_nodes = tada.taxdmp_merged_nodes(snakemake.input.merged_nodes)
df["taxid"] = df['taxid'].map(merged_nodes).fillna(df['taxid'])

# Find Taxids that is in eukaryota, viruses, other,
# and unclassified
for taxid in df["taxid"].unique():
    curr_taxid = taxid
    rank = ""
    while rank != "domain":
        rank = taxonomy[curr_taxid][1]
        curr_taxid = taxonomy[curr_taxid][0]
        if curr_taxid in [2759, 10239, 28384, 12908]:
            if snakemake.params.context == "required_genomes":
                print(f"{taxid} is not in Archaea or Bacteria, please change the file with required genomes")
            elif snakemake.params.context == "sampling_scheme":
                print(f"{taxid} in not in Archaea or Bacteria, please change the sampling scheme")
            else:
                print(f"{taxid} in not in Archaea or Bacteria")
            sys.exit(1)

# Export cleaned dataset
df.to_csv(snakemake.output.dataset, sep="\t", index=False)
