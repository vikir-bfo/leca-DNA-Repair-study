"""
HMM_transposing.py to transpose per-genome HMMER output into per-profile tables.

Input
-----
snakemake.input : an iterable of hmmsearch --tblout files, one per genome. The
    file stem is taken as the genome/species identifier, so file naming carries
    the sample metadata.

Output
------
One tab-separated, header-less file per profile, written to the directory of the
first declared Snakemake output, named <profile_id>.tsv with the columns
query, subject, species, evalue, score, bias, description.

"""

import os
import pandas as pd

input_dir = snakemake.input
first_output_folder = snakemake.output[0]
output_folder = os.path.dirname(first_output_file)
os.makedirs(output_folder, exist_ok=True)
cols = ['query', 'subject', 'species', 'evalue', 'score', 'bias', 'description']
all_data = []

for filepath in input_dir:
    species_name = os.path.splitext(os.path.basename(filepath))[0]
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split(None, 18)
            if len(parts) >= 19:
                all_data.append({
                    "query":       parts[0],
                    "subject":     parts[2],
                    "species":     species_name,
                    "evalue":      parts[4],
                    "score":       parts[5],
                    "bias":        parts[6],
                    "description": parts[18]
                })
df = pd.DataFrame(all_data)
df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')
df['score']  = pd.to_numeric(df['score'],  errors='coerce')
df = df.sort_values(by=["subject", "score"], ascending=[True, False])
#one output table per profile, header-less
for ko_id, group in df.groupby('subject'):
    filepath = os.path.join(output_folder, f"{ko_id}.tsv")
    group.to_csv(filepath, sep='\t', index=False, header=False)
