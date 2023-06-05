from ete3 import Tree
import pandas as pd
import bisect
import random
import sys


def calculate_distance(node):
    child_nodes = node.get_children()
    if len(child_nodes) > 1 and all([child.is_leaf() for child in child_nodes]):
        distance = sum([node.get_distance(child) for child in child_nodes])
        return (distance, node)

phylogeny_path = snakemake.input.phylogeny
metadata_path = snakemake.input.metadata
nr_taxa = snakemake.params.taxa
completeness = float(snakemake.params.completeness)
contamination = float(snakemake.params.contamination)
prune_method = snakemake.params.prune_method
random.seed(int(snakemake.params.seed))
output_metadata = snakemake.output.metadata
output_phylogeny = snakemake.output.phylogeny


# Read phylogeny
tree = Tree(phylogeny_path, format=1, quoted_node_names=True)

# Prune taxa from tree according to filters
df = pd.read_csv(metadata_path, sep="\t")
df[
    ["domain", "phylum", "class", "order", "family", "genus", "species"]
] = df.gtdb_taxonomy.str.split(";", expand=True)

df["domain"] = df["domain"].str.replace("d__", "")
df["phylum"] = df["phylum"].str.replace("p__", "")
df["class"] = df["class"].str.replace("c__", "")
df["order"] = df["order"].str.replace("o__", "")
df["family"] = df["family"].str.replace("f__", "")
df["genus"] = df["genus"].str.replace("g__", "")
df["species"] = df["species"].str.replace("s__", "")

df = df[df['accession'].isin(tree.get_leaf_names())]

df = df[
    (df.checkm_contamination <= contamination)
    & (df.checkm_completeness >= completeness)
]

clean_accessions = list(set(df['accession'].to_list()))

# If we ask for 0 taxa, we still require the output-files for downstream
# steps of the workflow.
if nr_taxa == 0:
    df = df[0:0]
    df.to_csv(output_metadata, sep="\t", index=False)
    with open(output_phylogeny, "w") as f:
        pass
    sys.exit()

orig_tree_size = len(set(tree.get_leaf_names()))
print(f"Tree size: {orig_tree_size}")
tree.prune(clean_accessions, preserve_branch_length=True)
filt_tree_size = len(set(tree.get_leaf_names()))
print(f"Tree size after filtering out taxa: {filt_tree_size}")

# Compute distance between each sister leaf-pair
print('Calculate pair-wise distances')
distance_list = []
nodes = tree.traverse()
for node in nodes:
    child_nodes = node.get_children()
    if len(child_nodes) > 1 and all([child.is_leaf() for child in child_nodes]):
        distance = sum([node.get_distance(child) for child in child_nodes])
        distance_list.append((distance, node))

# Sort the distances
sorted_distance_list = sorted(distance_list, key=lambda t: t[0])

# Until selected number of taxa
print('Start to remove leaves')
while len(tree.get_leaves()) > nr_taxa:
    if len(tree.get_leaves()) % 1000 == 0:
        print(f"Tree size  {len(tree.get_leaves())}")
    min_distance, min_node = sorted_distance_list[0]

    # Randomly trim one of the leaf
    min_node_children = min_node.get_children()
    if prune_method == "random":
        prune_idx = random.randint(0,1)
        if prune_idx == 0:
            keep_idx = 1
        else:
            keep_idx = 0
    elif prune_method == "shortest":
        if min_node_children[0].dist <= min_node_children[1].dist:
            keep_idx = 0
            prune_idx = 1
        else:
            keep_idx = 1
            prune_idx = 0
    elif prune_method == "longest":
        if min_node_children[0].dist >= min_node_children[1].dist:
            keep_idx = 0
            prune_idx = 1
        else:
            keep_idx = 1
            prune_idx = 0
    prune = min_node_children[prune_idx]
    keep = min_node_children[keep_idx]
    parent_dist = keep.dist
    new_parent_dist = parent_dist + keep.up.dist
    parent = keep.up
    new_parent = keep.up.up

    # Trim a leaf and remove the parent node
    prune.detach()
    if len(parent.get_children()) == 1:
        parent.delete()

    # Updated the distant to the parent node for the
    # node that remains.
    keep.dist = new_parent_dist

    # Remove the node from the list
    sorted_distance_list.pop(0)

    # Calculate new distance and insert into list
    child_nodes = new_parent.get_children()
    if len(child_nodes) > 1 and all([child.is_leaf() for child in child_nodes]):
        distance = sum([new_parent.get_distance(child) for child in child_nodes])
        bisect.insort_left(sorted_distance_list, (distance, new_parent), key=lambda i: i[0])

leafs = tree.get_leaf_names()

print("Write output files")
df = df[df["accession"].isin(leafs)]
df["accession"] = df["accession"].str.replace("GB_", "")
df["accession"] = df["accession"].str.replace("RS_", "")
df.to_csv(output_metadata, sep="\t", index=False)
tree.write(outfile=output_phylogeny)

