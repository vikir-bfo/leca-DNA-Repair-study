"""
Best_Unique_Hits.py is used to collapse the genome-level HMM hits to one best hit per gene.

Input
-----
Every *.tsv file in the working directory (excluding the output directory),
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

# Gene-grouping conventions, tried in order. Each entry is (match_first_token_only,
# The captured group becomes the gene-group key, so isoform and ORF suffixes
# are stripped and collapse to a single gene.
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

def gene_group(desc):
    """Return the gene-level grouping key from the hit description.
    """
    desc_str = str(desc).strip()
    if not desc_str:
        return None 
    first = desc_str.split()[0]
    for use_first, pat in _GROUP_PATTERNS:
        target = first if use_first else desc_str
        m = (re.match if use_first else re.search)(pat, target)
        if m:
            return m.group(1)
            
    # no gene information: recover the component identifier from elsewhere in the
    # description instead.
    if re.match(r'm\.\d+', first):
        for sub_pat in [r'(comp\d+_c\d+_seq\d+)', r'(Unigene\d+_\w+)']:
            m2 = re.search(sub_pat, desc_str)
            if m2:
                return m2.group(1)
    return None
                

def get_gene(desc):
    """Extract the annotated product name from a description, without any isoform
    qualifier.
    """
    match = re.search(r'^\S+\s+(.+?)(?=\s+isoform|$)', str(desc))
    if match:
        return match.group(1).strip()
    return str(desc).strip()

def process_file(input_files, output_path):
    """Collapse one or more hit tables to the best hit per gene and write the result.
    """
    combined_dfs = []
    for seq_path in input_files:
        if os.path.getsize(seq_path) == 0:
            print(f"  -> Skipping {seq_path} (File is empty)")
            continue
        try:
            df = pd.read_csv(seq_path, sep='\t', header=None, names=input_cols, comment='#')
            print(f"  -> Read {len(df)} rows from {seq_path.name}")
            combined_dfs.append(df)
        except pd.errors.EmptyDataError:
            print(f"  -> Skipping {seq_path} (No parseable data/only comments)")
            continue
        except Exception as e:
            print(f"  -> Error parsing {seq_path}: {e}")
            continue

    if not combined_dfs:
        print(f"  -> No valid data found for {output_path.name}. Skipping output.")
        return

    all_dfs = pd.concat(combined_dfs, ignore_index=True)
    print(f"  -> Total rows before processing: {len(all_dfs)}")
    
    all_dfs['description'] = all_dfs['description'].fillna('').astype(str)
    all_dfs['gene_group'] = all_dfs['description'].apply(gene_group)
    all_dfs['gene_group'] = all_dfs['gene_group'].fillna(all_dfs['seq_id'])
    all_dfs['gene_name'] = all_dfs['description'].apply(get_gene)
    all_dfs['desc_clean'] = all_dfs['description'].str.strip()
    all_dfs['gene_group'] = all_dfs['gene_group'].replace('', pd.NA).fillna(all_dfs['seq_id'])
    all_dfs['gene_name'] = all_dfs['gene_name'].replace('', pd.NA).fillna(all_dfs['seq_id'])
    all_dfs['desc_clean'] = all_dfs['desc_clean'].replace('', pd.NA).fillna(all_dfs['seq_id'])

    all_dfs['evalue'] = pd.to_numeric(all_dfs['evalue'], errors='coerce')
    all_dfs['score'] = pd.to_numeric(all_dfs['score'], errors='coerce')
    
    # Tie-breaking using the ORF log-likelihood ('score=') and the ORF length ('len:').
    # Absent values become 0 so that annotated ORFs outrank unannotated ones.
    all_dfs['orf_score'] = pd.to_numeric(
        all_dfs['description'].str.extract(r'score=([\d.]+)')[0], errors='coerce'
    ).fillna(0)
    
    all_dfs['orf_len'] = pd.to_numeric(
        all_dfs['description'].str.extract(r'len:(\d+)')[0], errors='coerce'
    ).fillna(0)
    
    # Step 4: E-value ascending first (primary criterion of homology support), then
    # bit score, ORF score and ORF length descending.
    all_sorted = all_dfs.sort_values(
        ['gene_group', 'evalue', 'score', 'orf_score', 'orf_len'],
        ascending=[True, True, False, False, False]
    )
    
    # Step 5: keep the top-ranked hit per gene group, then per gene name.
    best_df = all_sorted.drop_duplicates(subset=['gene_group'], keep='first')
    best_df = best_df.drop_duplicates(subset=['gene_name'], keep='first')
    
    filtered_hits = best_df.sort_values(by=['evalue', 'score'], ascending=[True, False])
    
    print(f"  -> Total rows after filtering: {len(filtered_hits)}")
    if len(filtered_hits) > 0:
        filtered_hits = filtered_hits[output_cols] 
        filtered_hits.to_csv(output_path, sep='\t', index=False)
        print(f"  -> Saved to {output_path}")
    else:
        print(f"  -> All rows were filtered out for {output_path.name}!")

# every *.tsv in the current directory is processed independently and
# written under best_HMMs/ with its original filename.
if __name__ == "__main__":
    input_dir = Path(".") 
    output_dir = Path("best_HMMs")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tsv_files = [f for f in input_dir.glob("*.tsv") if output_dir.name not in f.parts]
    
    for tsv_file in tsv_files:
        output_path = output_dir / tsv_file.name
        process_file([tsv_file], output_path)
