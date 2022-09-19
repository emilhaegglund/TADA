import pandas as pd
from Bio import SeqIO
import os
import gzip
import sys

taxa_path = sys.argv[1]
proteome_dir = sys.argv[2]
output = sys.argv[3]

taxid_accession_df = pd.read_csv(
    taxa_path, sep="\t", names=["assembly_accession", "taxid"]
)

data = []
for f in os.listdir(proteome_dir):
    f_path = os.path.join(proteome_dir, f)
    assembly_accession = "_".join(f.split("_")[:2])
    with gzip.open(f_path, "rt") as stream:
        for record in SeqIO.parse(stream, "fasta"):
            data.append([assembly_accession, record.id])

protein_df = pd.DataFrame(data, columns=["assembly_accession", "accession.version"])

df = pd.merge(taxid_accession_df, protein_df, on="assembly_accession")
df[["accession.version", "taxid"]].to_csv(output, sep='\t', index=False)
