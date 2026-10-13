"""
Best_Unique_Hits.py is used to collapse the genome-level HMM hits to one best hit per gene.

Input
-----
snakemake.input : List of *.tsv files in the working directory (excluding the output directory),
header-less, with the field order produced by HMM_transposing.py:
seq_id, query, subject, evalue, score, bias, description -- where seq_id is the
subject protein/transcript, query is the profile identifier and subject is the
genome or species label.

Output
------
best_HMMs/<same filename>, tab-separated with a header row and the columns
query, seq_id, evalue, score, bias, description.

"""
import re
import pandas as pd
import os
from pathlib import Path

# Field order of the header-less tables and of the written output.
# 'query' is the profile identifier; 'seq_id' is the hit protein/transcript.
input_cols = ['seq_id', 'query', 'subject', 'evalue', 'score', 'bias', 'description']
output_cols = ['query', 'seq_id', 'evalue', 'score', 'bias', 'description']
_GROUP_PATTERNS = [
    (True, r'(TRINITY_DN\d+_c\d+_g\d+)_i'),
    (True, r'(\w+_TRINITY_DN\d+_c\d+_g\d+)'),
    (True, r'(c\d+_g\d+)_i'),
    (True, r'(gb\|[^|]+\|)\.p\d+$'),
    (True, r'(.+?)\.p\d+$'),
    (True, r'(comp\d+_c\d+_seq\d+)\|'),
    (True, r'(Partr_v1_DN\d+_c\d+_g\d+)'),
    (False, r'(TR\d+_c\d+_g\d+)_i'),
]
# Each entry is (match_first_token_only,
# The captured group becomes the gene-group key, so isoform and ORF suffixes
# are stripped and collapse to a single gene.
def gene_group(desc):
    first = desc.split()[0]
    for use_first, pat in _GROUP_PATTERNS:
        target = first if use_first else desc
        m = (re.match if use_first else re.search)(pat, target)
        if m:
            return m.group(1)
            # no gene information: recover the component identifier in the
    # description instead.
    if re.match(r'm\.\d+', first):
        for sub_pat in [r'(comp\d+_c\d+_seq\d+)', r'(Unigene\d+_\w+)']:
            m2 = re.search(sub_pat, desc)
            if m2:
                return m2.group(1)
    return first

def get_gene(desc):
    """Extract the annotated product name from a description, without any isoform
    qualifier.
    """
    match = re.search(r'^\S+\s+(.+?)(?=\s+isoform|$)', str(desc))
    if match:
        return match.group(1).strip()
    return str(desc).strip()

def process_file(input_files, output_path):
    """Collapse one or more hit tables to the best hit per gene 
        and write the result.
    """
    combined_dfs = []
    for seq_path in input_files:
        df = pd.read_csv(seq_path, sep='\t', header=None, names=input_cols)
        combined_dfs.append(df)
    all_dfs = pd.concat(combined_dfs, ignore_index=True)
    all_dfs['gene_group'] = all_dfs['description'].apply(gene_group)
    all_dfs['gene_name'] = all_dfs['description'].apply(get_gene)
    all_dfs['desc_clean'] = all_dfs['description'].str.strip()
    all_dfs['evalue'] = pd.to_numeric(all_dfs['evalue'], errors='coerce')
    all_dfs['score'] = pd.to_numeric(all_dfs['score'], errors='coerce')
    all_dfs['orf_score'] = pd.to_numeric(
        all_dfs['description'].str.extract(r'score=([\d.]+)')[0], errors='coerce'
    ).fillna(0)
    # E-value ascending first (primary criterion of homology support), then
    # bit score, ORF score and ORF length descending.
    all_dfs['orf_len'] = pd.to_numeric(
        all_dfs['description'].str.extract(r'len:(\d+)')[0], errors='coerce'
    ).fillna(0)
    all_sorted = all_dfs.sort_values(
    ['gene_group', 'evalue', 'score', 'orf_score', 'orf_len'],
    ascending=[True, True, False, False, False])
    # keep the top-ranked hit per gene group, then per gene name.
    best_df = all_sorted.drop_duplicates(subset=['gene_group'], keep='first')
    best_df = best_df.drop_duplicates(subset=['gene_name'], keep='first')
    best_df = best_df.drop_duplicates(subset=['desc_clean'], keep='first')
    filtered_hits = best_df.sort_values(by=['evalue', 'score'], ascending=[True, False])
    filtered_hits.to_csv(output_path, sep='\t', index=False)
process_file(snakemake.input, snakemake.output[0])
