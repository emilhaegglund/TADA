import yaml
import pandas as pd
import sys


def check_taxa_name(taxa, taxa_levels, df):
    for level in taxa_levels:
        if taxa in df[level].to_list():
            return True
    return False


def get_taxa_level_index(taxa, taxa_levels, df):
    for taxa_level_index, level in enumerate(taxa_levels):
        if taxa in df[level].to_list():
            return taxa_level_index


with open("../../config/sampling_scheme_planctos.yaml", "r") as stream:
    sampling_scheme = yaml.safe_load(stream)
df = pd.read_csv("test.tsv", sep="\t", low_memory=False)
df["phylum"] = df["phylum"].str.replace("p__", "")
df["class"] = df["class"].str.replace("c__", "")
df["order"] = df["order"].str.replace("o__", "")
df["family"] = df["family"].str.replace("f__", "")
df["genus"] = df["genus"].str.replace("g__", "")
df["species"] = df["species"].str.replace("s__", "")

taxa_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
sampling_order = {}
for taxa in sampling_scheme:
    if check_taxa_name(taxa, taxa_levels, df):
        sampling_level = sampling_scheme[taxa]["sampling_level"]
        n_taxa = sampling_scheme[taxa]["taxa"]
        taxa_level_index = get_taxa_level_index(taxa, taxa_levels, df)
        sampling_level_index = taxa_levels.index(sampling_level)
        # Make sure that we not sample from a higher taxonomic level
        # compared to the given taxonomic name
        if sampling_level_index >= taxa_level_index:
            if taxa_level_index in sampling_order.keys():
                sampling_order[taxa_level_index].append(
                    [taxa, sampling_level_index, n_taxa]
                )
            else:
                sampling_order[taxa_level_index] = [[taxa, sampling_level, n_taxa]]
        else:
            sys.exit("Can't sample from a higher taxonomic level")
    else:
        sys.exit("{taxa} is not present in GTDB")

# Reorder the sampling dictionary so tha we start sampling from species level and the continue
# with higher levels
sampling_order = {
    key: sampling_order[key] for key in sorted(sampling_order.keys(), reverse=True)
}
sampled_dfs = []  # Store sampled records
used_data = []
pd.set_option('display.max_colwidth', None)
for taxa_level_index in sampling_order.keys():
    taxa_level = taxa_levels[taxa_level_index]
    print(taxa_level)
    for sampling in sampling_order[taxa_level_index]:
        taxa = sampling[0]
        sampling_level = sampling[1]
        print(sampling_level)
        n_taxa = sampling[2]
        sampling_df = df[df[taxa_level] == taxa]
        for used_df in used_data:
            sampling_df = pd.concat([sampling_df, used_df]).drop_duplicates(
                keep=False
            )
        for i, taxa_level_df in sampling_df.groupby(sampling_level):
            print(taxa_level_df[["taxonomy"]])
            # Remove data that has already been used to sample from.

            if taxa_level_df.shape[0] > n_taxa:
                sampled_dfs.append(taxa_level_df.sample(n_taxa))
            else:
                sampled_dfs.append(taxa_level_df)
        used_data.append(sampling_df)

sampled_df = pd.concat(sampled_dfs)
sampled_df[["assembly_accession_x", "taxonomy"]].to_csv(
    "sampling_scheme_planctos.output.tsv", sep="\t", index=False
)
# for i, taxa_level_df in df.groupby(args.taxonomic_level):
# Can't take a sample if the sample size we ask for is larger than
# the number of taxa in that group. In that case, use all taxa in
# the group.
#    if taxa_level_df.shape[0] > args.max_taxa:
#        sampled_accessions.append(taxa_level_df.sample(args.max_taxa))
#    else:
#        sampled_accessions.append(taxa_level_df)
