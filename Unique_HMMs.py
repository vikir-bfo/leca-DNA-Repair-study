"""
Unique_HMMs.py to reduce HMM hit tables to one best hit per subject sequence and
remove hits belonging to the Asgard lineage.
"""

import pandas as pd
import os
from pathlib import Path

output_cols = ['query', 'seq_id', 'evalue', 'score', 'bias', 'description']

PROTTOTAXID_PATH = "/local/two/vineeth/nr_cluster/ncbi_taxonomyasgard_prot.accession2taxid.tsv"
LINEAGE_PATH = "/local/two/vineeth/DNArep_asCOG-nr/gene_trees/focal_subTree/fullnamelineage.dmp"

def load_prot_to_taxid(prottotaxid_path):
    """Map protein accession to NCBI taxon identifier.
    """
    prot_to_taxid = {}
    try:
        with open(prottotaxid_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    accession = parts[0]
                    taxid = parts[2]
                    prot_to_taxid[accession] = taxid
    except Exception as e:
        print(f"  -> Warning: Could not load prottotaxid file: {e}")
    return prot_to_taxid

def load_lineage_dict(lineage_path):
    """Map NCBI taxon identifier to its full lineage string.
    """
    lineage_dict = {}
    try:
        with open(lineage_path, 'r') as f:
            for line in f:
                parts = [p.strip() for p in line.split('|')]
                if len(parts) >= 3:
                    tax_id = parts[0]
                    lineage = parts[2]
                    lineage_dict[tax_id] = lineage
    except Exception as e:
        print(f"  -> Warning: Could not load lineage file: {e}")
    return lineage_dict

def is_prometharchaea(protein_id, prot_to_taxid, lineage_dict):
    """Test whether a subject accession belongs to the Asgard lineage.
    """
    if not prot_to_taxid or not lineage_dict:
        return False
    if protein_id not in prot_to_taxid:
        return False
    taxid = prot_to_taxid[protein_id]
    if taxid not in lineage_dict:
        return False
    lineage = lineage_dict[taxid]
    return 'Prometharchaeota' in lineage

def process_file(input_files, output_path, prot_to_taxid=None, lineage_dict=None):
    """Filter and de-duplicate the hit tables, writing the best hit per
    subject sequence.
    """
    parsed_data = []
    prometharchaea_count = 0
    
    for seq_path in input_files:
        if os.path.getsize(seq_path) == 0:
            print(f"  -> Skipping {seq_path} (File is empty)")
            continue
        try:
            with open(seq_path, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split()
                    if len(parts) < 17:
                        continue
                    
                    seq_id = parts[1]
                    
                    if is_prometharchaea(seq_id, prot_to_taxid, lineage_dict):
                        prometharchaea_count += 1
                        continue
                    
                    parsed_data.append({
                        'query': parts[2],
                        'seq_id': seq_id, 
                        'evalue': parts[3],
                        'score': parts[4],
                        'bias': parts[5],
                        'description': " ".join(parts[17:])
                    })
        except Exception as e:
            print(f"  -> Error parsing {seq_path}: {e}")
            continue

    all_dfs = pd.DataFrame(parsed_data)
    all_dfs['evalue'] = pd.to_numeric(all_dfs['evalue'], errors='coerce')
    all_dfs['score'] = pd.to_numeric(all_dfs['score'], errors='coerce')
    all_sorted = all_dfs.sort_values(by=['evalue', 'score'], ascending=[True, False])
    best_df = all_sorted.drop_duplicates(subset=['seq_id'], keep='first')
    if len(best_df) > 0:
        filtered_hits = best_df[output_cols] 
        filtered_hits.to_csv(output_path, sep='\t', index=False)
    else:
        print(f"  -> All rows were filtered out for {output_path.name}!")

if __name__ == "__main__":
    prot_to_taxid = load_prot_to_taxid(PROTTOTAXID_PATH)
    lineage_dict = load_lineage_dict(LINEAGE_PATH)
    
    input_dir = Path(".") 
    output_dir = Path("best_HMMs")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tsv_files = [f for f in input_dir.glob("*.tsv") if output_dir.name not in f.parts]
    
    for tsv_file in tsv_files:
        output_path = output_dir / tsv_file.name
        process_file([tsv_file], output_path, prot_to_taxid, lineage_dict)
