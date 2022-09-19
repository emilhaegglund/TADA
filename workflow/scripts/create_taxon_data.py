import pandas as pd
import argparse


def parse_command_line():
    args = argparse.ArgumentParser()
    args.add_argument("--metadata", required=True)
    args.add_argument("--out-nodes", required=True)
    args.add_argument("--out-names", required=True)
    args.add_argument("--out-taxid", required=True)

    return args.parse_args()


args = parse_command_line()

df = pd.read_csv(
    args.metadata,
    sep="\t",
)

taxa_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# First, drop the columns with out rank-prefix
# This has to be done because some taxa in different ranks
# have the same name, for example there is both an order called
# UBA2392 and a family called UBA2392.
df.drop(axis="columns", labels=taxa_levels, inplace=True)

# Then split the taxonomy column to new columns
df[
    ["domain", "phylum", "class", "order", "family", "genus", "species"]
] = df.gtdb_taxonomy.str.split(";", expand=True)

# Start taxids from 2, 1 is reserved for root
counter = 2
taxid_name_map = {}  # {name:taxid}
names_out = open(args.out_names, "w")

# Give all taxa at all ranks a taxid
for taxa_level in taxa_levels:
    taxas = set(df[taxa_level].to_list())
    for taxa in taxas:
        taxid_name_map[taxa] = counter
        names_out.write(
            str(counter) + "\t|\t" + taxa + "\t|\t" + "" + "\t|\t" + "GTDB\t|\n"
        )  # write taxids to names file
        counter += 1

# Give all accessions a taxid
for accession in set(df["accession"].to_list()):
    names_out.write(
        str(counter) + "\t|\t" + accession + "\t|\t" + "" + "\t|\t" + "GTDB\t|\n"
    )  # write taxids to names file
    taxid_name_map[accession] = counter
    counter += 1

# Write an accession to taxid map
with open(args.out_taxid, "w") as f:
    for accession in set(df["accession"].to_list()):
        f.write(accession + "\t" + str(taxid_name_map[accession]) + "\n")
names_out.close()

# Build a dictionary which contain the taxid as key and the parent node taxid
# and the rank as values for each node. The root-node is already given
parent_taxid_map = {1: [1, "no rank"]}
for i, row in df.iterrows():
    for level in ["accession"] + taxa_levels[::-1]:
        if level == "accession":
            prev_taxa = row[level]
            prev_rank = "isolate"
        else:
            curr_taxa = row[level]
            curr_rank = level
            parent_taxid = taxid_name_map[curr_taxa]
            taxid = taxid_name_map[prev_taxa]
            if taxid not in parent_taxid_map.keys():
                parent_taxid_map[taxid] = [parent_taxid, prev_rank]
            if curr_rank == "domain":
                if 2 not in parent_taxid_map.keys():
                    parent_taxid_map[2] = [1, "superkingdom"]
            prev_taxa = curr_taxa
            prev_rank = curr_rank

# Write the nodes to file. The output has to be sorted according to
# the taxids.
with open(args.out_nodes, "w") as f:
    for node in sorted(parent_taxid_map.keys()):
        if node == 1:  # Root is special case
            node_str = "\t|\t".join(
                [
                    str(node),
                    str(parent_taxid_map[node][0]),
                    parent_taxid_map[node][1],
                    "",
                    "8",
                    "0",
                    "1",
                    "0",
                    "0",
                    "0",
                    "0",
                    "0",
                ]
            )
        else:
            node_str = "\t|\t".join(
                [
                    str(node),
                    str(parent_taxid_map[node][0]),
                    parent_taxid_map[node][1],
                    "XX",
                    "0",
                    "1",
                    "11",
                    "1",
                    "0",
                    "1",
                    "1",
                    "0",
                ]
            )
        f.write(node_str + "\t|\n")
