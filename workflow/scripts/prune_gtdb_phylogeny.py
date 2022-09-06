from ete3 import Tree
import argparse
import pandas as pd
import time
import bisect


def parse_command_line():

    args = argparse.ArgumentParser()
    args.add_argument("--phylogeny", required=True)
    args.add_argument("--gtdb-metadata", required=True)
    args.add_argument("--taxa", required=True, type=int)
    args.add_argument("--output-refseq", required=True)
    args.add_argument("--output-genbank", required=True)
    args.add_argument("--completeness", type=float, default=0)
    args.add_argument("--contamination", type=float, default=100)
    args.add_argument("--cpus", default=1, type=int)
    return args.parse_args()


args = parse_command_line()

def calculate_distance(node):
    child_nodes = node.get_children()
    if len(child_nodes) > 1 and all([child.is_leaf() for child in child_nodes]):
        distance = sum([node.get_distance(child) for child in child_nodes])
        return (distance, node)

# Read phylogeny
tree = Tree(args.phylogeny, format=1, quoted_node_names=True)



# First prune taxa from tree according to filter
df = pd.read_csv(args.gtdb_metadata, sep="\t")

df = df[df['accession'].isin(tree.get_leaf_names())]
print(df.shape)

df = df[
    (df.checkm_contamination <= args.contamination)
    & (df.checkm_completeness >= args.completeness)
]

clean_accessions = list(set(df['accession'].to_list()))
print(len(clean_accessions))

start_time = time.time()
print(len(set(tree.get_leaf_names())))
tree.prune(clean_accessions, preserve_branch_length=True)
print(len(tree.get_leaf_names()))
print(time.time() - start_time)

# Compute distance between each sister leaf-pair
print('Calculate pair-wise distances')
distance_list = []
nodes = tree.traverse()
for node in nodes:
    child_nodes = node.get_children()
    if len(child_nodes) > 1 and all([child.is_leaf() for child in child_nodes]):
        distance = sum([node.get_distance(child) for child in child_nodes])
        distance_list.append((distance, node))

#start_time = time.time()
#pool = mp.Pool(args.cpus)
#distance_list = pool.map(calculate_distance, nodes)
distance_list = [dist for dist in distance_list if dist != None]
#print(time.time() - start_time)

# Sort the distances
sorted_distance_list = sorted(distance_list, key=lambda t: t[0])
print(sorted_distance_list[:10])

# Until selected number of taxa
print('Start to remove leaves')
while len(tree.get_leaves()) >= args.taxa:
    if len(tree.get_leaves()) % 1000 == 0:
        print(len(tree.get_leaves()))
    min_distance, min_node = sorted_distance_list[0]

    # Randomly trim one of the leaf
    min_node_children = min_node.get_children()
    prune = min_node_children[1]
    keep = min_node_children[0]
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
        #sorted_distance_list.append((distance, new_parent))
        #sorted_distance_list = sorted(sorted_distance_list, key=lambda t: t[0])

leafs = tree.get_leaf_names()


df = df[df["accession"].isin(leafs)]
df["accession"] = df["accession"].str.replace("GB_", "")
df["accession"] = df["accession"].str.replace("RS_", "")
sampled_df_genbank = df[df["accession"].str.contains("GCA")]
sampled_df_refseq = df[df["accession"].str.contains("GCF")]
sampled_df_refseq["accession"].to_csv(args.output_refseq, sep="\t", index=False)
sampled_df_genbank["accession"].to_csv(args.output_genbank, sep="\t", index=False)

