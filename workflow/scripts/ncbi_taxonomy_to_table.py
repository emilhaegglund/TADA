from ete3 import NCBITaxa
import pandas as pd
import argparse


def parse_command_line():

    args = argparse.ArgumentParser()
    args.add_argument("--assembly-summary", required=True)
    args.add_argument("--names", required=True)
    args.add_argument("--nodes", required=True)
    args.add_argument("--taxdmp", required=True)
    args.add_argument("--output", required=True)

    return args.parse_args()


args = parse_command_line()

# Use the NCBITaxa from ete to find the lineage of each taxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database(args.taxdmp)

header = [
    "assembly_accession",
    "bioproject",
    "biosample",
    "wgs_master",
    "refseq_category",
    "taxid",
    "species_taxid",
    "organism_name",
    "infraspecific_name",
    "isolate",
    "version_status",
    "assembly_level",
    "release_type",
    "genome_rep",
    "seq_rel_date",
    "asm_name",
    "submitter",
    "gbrs_paired_asm",
    "paired_asm_comp",
    "ftp_path",
    "excluded_from_refseq",
    "relation_to_type_material",
    "asm_not_live_date",
]

assemblies_df = pd.read_csv(args.assembly_summary, sep="\t", comment="#", names=header)

# Create a dataframe from the names.dmp file
names_dict = {"taxid": [], "name_txt": [], "uniq_name": [], "name_class": []}
with open(args.names, "r") as f:
    for line in f:
        line = line.strip()
        line = line.split("|")
        line = [l.strip() for l in line]
        names_dict["taxid"].append(line[0])
        names_dict["name_txt"].append(line[1])
        names_dict["uniq_name"].append(line[2])
        names_dict["name_class"].append(line[3])

names_df = pd.DataFrame.from_dict(names_dict)

# Create a dataframe from the nodes.dmp file
taxa_dict = {"taxid": [], "rank": []}
with open(args.nodes, "r") as f:
    for line in f:
        line = line.strip()
        line = line.split("|")
        line = [l.strip() for l in line]
        taxa_dict["taxid"].append(line[0]),
        taxa_dict["rank"].append(line[2])

taxa_df = pd.DataFrame.from_dict(taxa_dict)

# Combine names and taxa data frames
names_df = pd.merge(left=names_df, right=taxa_df, left_on="taxid", right_on="taxid")

# Use only the standard ranks
std_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
names_df = names_df[names_df["rank"].isin(std_ranks)]

# Store taxonomy data for each species
species_taxonomy = []
for taxid in assemblies_df["species_taxid"].unique():
    lineage = ncbi.get_lineage(taxid)
    lineage = [str(l) for l in lineage]
    lineage_name_df = names_df[names_df["taxid"].isin(lineage)]
    assemblies = assemblies_df[assemblies_df["species_taxid"] == taxid]
    names = lineage_name_df[lineage_name_df["name_class"] == "scientific name"][
        ["name_txt", "rank"]
    ]
    names = names.T
    names.columns = names.iloc[1]
    names = names.drop(names.index[1])
    for rank in std_ranks:
        if rank in names.columns:
            assemblies[rank] = names[rank].values[0]
        # If standard rank is missing, add NA
        else:
            assemblies[rank] = "NA"
    species_taxonomy.append(assemblies[["assembly_accession"] + std_ranks])

species_taxonomy = pd.concat(species_taxonomy)
species_taxonomy.to_csv(args.output, sep="\t", index=False)
