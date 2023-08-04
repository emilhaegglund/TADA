import pandas as pd
import sys


def read_taxonomy(nodes_path):
    taxonomy = {}
    with open(nodes_path, "r") as f:
        for line in f:
            line = line.strip("\n")
            line = line.split("|")
            taxid = int(line[0].strip("\t"))
            parent_taxid = int(line[1].strip("\t"))
            rank = line[2].strip("\t")
            if rank == "superkingdom":
                rank = "domain"
            taxonomy[taxid] = [parent_taxid, rank]

    return taxonomy

# Read original dataset
df = pd.read_csv(snakemake.input.dataset, sep="\t")
if "Organism Taxonomic ID" in df.columns:
    df.rename(columns={"Organism Taxonomic ID":"taxid"}, inplace=True)

# Read NCBI taxonomy
taxonomy = read_taxonomy(snakemake.input.nodes)

# Find Taxids that is in eukaryota, viruses, other,
# and unclassified
taxid_to_exclude = []
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
