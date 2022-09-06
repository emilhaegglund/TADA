import sys
import os


all_accessions_file = sys.argv[1]
proteome_dir = sys.argv[2]
output = sys.argv[3]

all_accessions = []
with open(all_accessions_file, 'r') as fin:
    for line in fin:
        line = line.strip()
        all_accessions.append(line)

downloaded_accessions = []
for f in os.listdir(proteome_dir):
    accession = "_".join(f.split("_")[:2])
    downloaded_accessions.append(accession)

with open(output, 'w') as fout:
    for accession in all_accessions:
        if accession not in downloaded_accessions:
            fout.write(accession + '\n')
