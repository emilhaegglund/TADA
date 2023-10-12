from types import NoneType
import pandas as pd
import yaml
import tada
import sys

# Define variables from snakemake
gtdb_metadata = snakemake.input.metadata
required_genomes_path = snakemake.input.required_genomes
sampling_scheme_path = snakemake.params.sampling_scheme
completeness = float(snakemake.params.completeness)
contamination = float(snakemake.params.contamination)
gtdb_representative = snakemake.params.gtdb_representative
seed = int(snakemake.params.seed)
output = snakemake.output[0]

# Read sampling scheme
with open(sampling_scheme_path, "r") as stream:
    sampling_scheme = yaml.safe_load(stream)

# Read GTDB metadata
df = pd.read_csv(gtdb_metadata, sep="\t", low_memory=False)
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


# Remove accession prefix
df["accession"] = df["accession"].str.replace("GB_", "")
df["accession"] = df["accession"].str.replace("RS_", "")

# Extract data for required genomes, this has to be done before all
# other filtering steps.
# If no required file has been given, the input type is not string
if type(required_genomes_path) == str:
    required_genomes_df = pd.read_csv(required_genomes_path, sep="\t")
    required_accessions = required_genomes_df["assembly_accession"].to_list()
    all_accessions = df["accession"].to_list()
    for accession in required_accessions:
        if accession not in all_accessions:
            print(f"{accession} not in GTDB")
            sys.exit(1)
    required_genomes_df = pd.merge(
        left=required_genomes_df,
        right=df,
        left_on="assembly_accession",
        right_on="accession",
    )


# Filter based on contamination and completeness
df = df[
    (df.checkm_contamination <= contamination)
    & (df.checkm_completeness >= completeness)
]

# Keep only GTDB-species representatives
if gtdb_representative:
    df = df[df.gtdb_representative == "t"]


# Define taxa levels
taxa_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
sampling_order = {}  # Store sampling parameters

if sampling_scheme is None:
    sampling_scheme = {}
if "all" in sampling_scheme.keys():
    sampling_parameters = sampling_scheme["all"]
    del sampling_scheme["all"]
    sampling_scheme["Bacteria"] = sampling_parameters
    sampling_scheme["Archaea"] = sampling_parameters

for taxa in sampling_scheme:
    if tada.check_taxa_name(taxa, taxa_levels, df):
        sampling_level = sampling_scheme[taxa]["sampling_level"]
        n_taxa = sampling_scheme[taxa]["taxa"]
        taxa_level_index = tada.get_taxa_level_index(taxa, taxa_levels, df)
        sampling_level_index = taxa_levels.index(sampling_level)
        # Make sure that we not sample from a higher taxonomic level
        # compared to the given taxonomic name
        if sampling_level_index >= taxa_level_index:
            # Store sampling parameters to dictionary
            if taxa_level_index in sampling_order.keys():
                sampling_order[taxa_level_index].append([taxa, sampling_level, n_taxa])
            else:
                sampling_order[taxa_level_index] = [[taxa, sampling_level, n_taxa]]
        else:
            sys.exit("Can't sample from a higher taxonomic level")
    else:
        sys.exit(f"{taxa} is not present in GTDB")

# Reorder the sampling dictionary so that we start sampling from species level and the continue
# with higher levels
sampling_order = {
    key: sampling_order[key] for key in sorted(sampling_order.keys(), reverse=True)
}
sampled_dfs = []  # Store sampled records
used_data = []  # Store data that has already been used to sample from.
for taxa_level_index in sampling_order.keys():
    taxa_level = taxa_levels[taxa_level_index]
    for sampling in sampling_order[taxa_level_index]:
        # Extract sampling parameters
        taxa = sampling[0]
        sampling_level = sampling[1]
        n_taxa = sampling[2]
        if n_taxa == "all":
            n_taxa = 0

        # Create a dataframe to sample from based on the selected taxa
        sampling_df = df[df[taxa_level] == taxa]

        # Remove data that has already been used to sample fromh
        exclude_entries = []
        for used_df in used_data:
            exclude_entries += used_df["accession"].to_list()
        sampling_df = sampling_df[~sampling_df["accession"].isin(exclude_entries)]

        # Find level under sampling level
        if sampling_level != "species":
            index = taxa_levels.index(sampling_level)
            index += 1
            base_level = taxa_levels[index]
            n_units = len(set(sampling_df[base_level].to_list()))
            prob_units = 1 / n_units
            dfs = []
            # Assign sampling probabilities to each taxa
            for i, taxa_level_df in sampling_df.groupby(base_level):
                n_base_level_taxa = taxa_level_df.shape[0]
                prob_base_level_taxa = prob_units / n_base_level_taxa
                taxa_level_df["sampling_prob"] = prob_base_level_taxa
                dfs.append(taxa_level_df)
            sampling_df = pd.concat(dfs)

        # Group the sampling dataframe based on the sampling level
        for i, taxa_level_df in sampling_df.groupby(sampling_level):
            sample_n_taxa = n_taxa
            # First extract and count required genomes
            if type(snakemake.input.required_genomes) == str:
                tmp_df = required_genomes_df[required_genomes_df[sampling_level] == i]
                if tmp_df.shape[0] > 0:

                    # Update the number of genomes that must be sampled
                    # after required genomes have been included and add
                    # required genomes to the sampled dataframe.
                    sampled_dfs.append(tmp_df)
                    sample_n_taxa = n_taxa - tmp_df.shape[0]

                    # Remove required genomes from taxa_level_df so that
                    # they are not included in the sampling.
                    tmp_df_accessions = tmp_df["accession"].to_list()
                    taxa_level_df = taxa_level_df[
                        ~taxa_level_df["accession"].isin(tmp_df_accessions)
                    ]
                    # Skip sampling if the number of required genomes are covering
                    # the number of genomes that should have been sampled.
                    if sample_n_taxa <= 0:
                        continue

            # Can't take a sample if the sample size we ask for is larger than
            # the number of taxa in that group. In that case, use all taxa in
            # the group.
            if sample_n_taxa == 0:
                sampled_dfs.append(taxa_level_df)
            elif taxa_level_df.shape[0] > sample_n_taxa:
                if "sampling_prob" in taxa_level_df.columns:
                    sampled_df = taxa_level_df.sample(
                        n=sample_n_taxa, weights="sampling_prob", random_state=seed
                    )
                    sampled_dfs.append(sampled_df)
                else:
                    sampled_df = taxa_level_df.sample(sample_n_taxa, random_state=seed)
                    sampled_dfs.append(sampled_df)
            else:
                sampled_dfs.append(taxa_level_df)
        # Finally, add the sampling dataframe to used data
        used_data.append(sampling_df)

# Also add the required genomes that are outside of the sampling scheme
if type(snakemake.input.required_genomes) == str:
    sampled_dfs.append(required_genomes_df)

sampled_df = pd.concat(sampled_dfs)
sampled_df.drop_duplicates(inplace=True)
sampled_df.to_csv(output, sep="\t", index=False)
