"""
tada

Description
-----------
This modules provides functions specific to the TADA workflow.

Functions
---------
taxdmp_names(names_path)
    Description: Read the 'names.dmp' file from NCBI-taxdmp and
        extract scientific names.
    Arguments:
        - names_path (str): Path to the 'names.dmp' file.
    Returns:
        - names (dict): A dictionary where keys are taxids and
            values are scientific names.

taxdmp_merged_nodes(merged_path)
    Description: Read the 'merged.dmp' file from NCBI-taxdmp and
        extract information about merged taxonomic nodes.
    Arguments:
        - merged_path (str): Path to the 'merged.dmp' file.
    Returns:
        - merged_nodes (dict): A dictionary where keys are the original
            taxids and the values are the merged taxids.

taxdmp_taxonomy(nodes_path)
    Description: Read the 'nodes.dmp' file from NCBI-taxdmp and
        extract information about the taxonomic hierarchy, including
        parent taxids and ranks.
    Arguments:
        - nodes_path (str): Path to the 'nodes.dmp' file.
    Returns:
        - taxonomy (dict): A dictionary where the keys are taxids, and
            values are lists containing parent taxid and the taxonomic
            rank.
"""


def taxdmp_names(names_path):
    """Read the 'names.dmp' file from NCBI-taxdmp and extract scientific names.

    Arguments:
    names_path: Path to names.dmp file

    Returns:
    names: Dictionary on form {taxid:scientific name}
    """
    names = {}
    with open(names_path, "r") as f:
        for line in f:
            line = line.strip("\n")
            line = line.split("\t|\t")
            if line[-1].strip("\t|") == "scientific name":
                names[int(line[0])] = line[1]
    return names


def taxdmp_merged_nodes(merged_path):
    """Read the merged.dmp file from NCBI-taxdmp

    Arguments:
    merged_path: Path to merged.dmp file

    Returns:
    merged_nodes: Dictionary on the form {taxid:merged_taxid}
    """
    merged_nodes = {}
    with open(merged_path, "r") as f:
        for line in f:
            line = line.strip("\n")
            line = line.strip("\t|")
            line = line.split("\t|\t")
            merged_nodes[int(line[0])] = line[1]

    return merged_nodes


def taxdmp_taxonomy(nodes_path):
    """Read the nodes.dmp file from NCBI-taxdmp and
    returns a dictionary with a list containing the
    taxid of the parent node and the taxonomic rank
    of the node.

    Arguments:
    nodes_path: Path to nodes.dmp file

    Return:
    taxonomy: Dictionary on the form {taxid:[parent_taxid, rank]}

    """
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


def check_taxa_name(taxa, taxa_levels, df):
    """
    Check that taxa is a valid taxa name
    """
    for level in taxa_levels:
        if taxa in df[level].to_list():
            return True
    return False


def get_taxa_level_index(taxa, taxa_levels, df):
    """
    Find the taxa level index for a taxa
    """
    for taxa_level_index, level in enumerate(taxa_levels):
        if taxa in df[level].to_list():
            return taxa_level_index
