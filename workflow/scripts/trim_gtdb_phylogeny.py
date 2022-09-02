from ete3 import Tree
import sys
import argparse


def parse_command_line():

    args = argparse.ArgumentParser()
    args.add_argument("--phylogeny", required=True)
    args.add_argument("--taxa", required=True, type=int)
    args.add_argument("--output", required=True)

    return args.parse_args()


args = parse_command_line()

# Read phylogeny
tree = Tree(args.phylogeny, format=1, quoted_node_names=True)

# Compute distance between each sister leaf-pair
print('Calculate pair-wise distances')
distance_list = []
for node in tree.traverse():
    child_nodes = node.get_children()
    if len(child_nodes) > 1 and all([child.is_leaf() for child in child_nodes]):
        distance = sum([node.get_distance(child) for child in child_nodes])
        distance_list.append((distance, node))

# Sort the distances
sorted_distance_list = sorted(distance_list, key=lambda t: t[0])

# Until selected number of taxa
print('Start to remove leaves')
while len(tree.get_leaves()) >= args.taxa:
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
        sorted_distance_list.append((distance, new_parent))
        sorted_distance_list = sorted(sorted_distance_list, key=lambda t: t[0])

tree.write(outfile=args.output)
