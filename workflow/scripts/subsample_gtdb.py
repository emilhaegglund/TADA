import pandas as pd
import argparse
import yaml
import sys


def read_command_line():
    """Read command line arguments"""
    args = argparse.ArgumentParser()
    args.add_argument("--gtdb-taxonomy", required=True)
    args.add_argument("--gtdb-metadata", required=True)
    args.add_argument("--sampling-scheme", required=True)
    args.add_argument("--output", required=True)
    args.add_argument("--completeness", type=float, default=0)
    args.add_argument("--contamination", type=float, default=100)
    args.add_argument("--gtdb-representative")

    return args.parse_args()


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

args = read_command_line()

# Read sampling scheme
with open(args.sampling_scheme, "r") as stream:
    sampling_scheme = yaml.safe_load(stream)

# Read GTDB taxonomy
taxa_df = pd.read_csv(
    args.gtdb_taxonomy, sep="\t", names=["assembly_accession", "taxonomy"]
)
taxa_df[
    ["domain", "phylum", "class", "order", "family", "genus", "species"]
] = taxa_df.taxonomy.str.split(";", expand=True)

taxa_df["domain"] = taxa_df["domain"].str.replace("d__", "")
taxa_df["phylum"] = taxa_df["phylum"].str.replace("p__", "")
taxa_df["class"] = taxa_df["class"].str.replace("c__", "")
taxa_df["order"] = taxa_df["order"].str.replace("o__", "")
taxa_df["family"] = taxa_df["family"].str.replace("f__", "")
taxa_df["genus"] = taxa_df["genus"].str.replace("g__", "")
taxa_df["species"] = taxa_df["species"].str.replace("s__", "")

# Read GTDB metadata
metadata_df = pd.read_csv(args.gtdb_metadata, sep="\t", low_memory=False)

# Merge the tables
df = pd.merge(
    left=taxa_df,
    right=metadata_df,
    left_on="assembly_accession",
    right_on="accession",
)

print(df.shape)
# Filter based on contamination and completeness
df = df[
    (df.checkm_contamination <= args.contamination)
    & (df.checkm_completeness >= args.completeness)
]

if args.gtdb_representative == True:
    df = df[df.gtdb_representative == "t"]


# Remove accession prefix
df["accession"] = df["accession"].str.replace("GB_", "")
df["accession"] = df["accession"].str.replace("RS_", "")

# Define taxa levels
taxa_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
sampling_order = {}  # Store sampling parameters

if 'all' in sampling_scheme.keys():
    sampling_parameters = sampling_scheme["all"]
    del(sampling_scheme['all'])
    sampling_scheme['Bacteria'] = sampling_parameters
    sampling_scheme['Archaea'] = sampling_parameters

for taxa in sampling_scheme:
    if check_taxa_name(taxa, taxa_levels, df):
        sampling_level = sampling_scheme[taxa]["sampling_level"]
        n_taxa = sampling_scheme[taxa]["taxa"]
        taxa_level_index = get_taxa_level_index(taxa, taxa_levels, df)
        sampling_level_index = taxa_levels.index(sampling_level)
        # Make sure that we not sample from a higher taxonomic level
        # compared to the given taxonomic name
        if sampling_level_index >= taxa_level_index:
            # Store sampling parameters to dictionary
            if taxa_level_index in sampling_order.keys():
                sampling_order[taxa_level_index].append(
                    [taxa, sampling_level, n_taxa]
                )
            else:
                sampling_order[taxa_level_index] = [[taxa, sampling_level, n_taxa]]
        else:
            sys.exit("Can't sample from a higher taxonomic level")
    else:
        print(taxa)
        sys.exit("{taxa} is not present in GTDB")

# Reorder the sampling dictionary so tha we start sampling from species level and the continue
# with higher levels
sampling_order = {
    key: sampling_order[key] for key in sorted(sampling_order.keys(), reverse=True)
}
sampled_dfs = []  # Store sampled records
used_data = []  # Store data that has already been used to sample from.
print(sampling_order)
for taxa_level_index in sampling_order.keys():
    taxa_level = taxa_levels[taxa_level_index]
    for sampling in sampling_order[taxa_level_index]:
        # Extract sampling parameters
        taxa = sampling[0]
        sampling_level = sampling[1]
        n_taxa = sampling[2]

        # Create a dataframe to sample from based on the selected taxa
        sampling_df = df[df[taxa_level] == taxa]

        # Remove data that has already been used to sample fromh
        for used_df in used_data:
            sampling_df = pd.concat([sampling_df, used_df]).drop_duplicates(keep=False)

        # Group the sampling dataframe based on the sampling level
        for i, taxa_level_df in sampling_df.groupby(sampling_level):

            # Can't take a sample if the sample size we ask for is larger than
            # the number of taxa in that group. In that case, use all taxa in
            # the group.
            if taxa_level_df.shape[0] > n_taxa:
                sampled_dfs.append(taxa_level_df.sample(n_taxa))
            else:
                sampled_dfs.append(taxa_level_df)

        # Finally, add the sampling dataframe to used data
        used_data.append(sampling_df)

sampled_df = pd.concat(sampled_dfs)
sampled_df.to_csv(args.output, sep='\t', index=False)
