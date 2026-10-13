from Bio import SeqIO
import pickle
from ete4 import PhyloTree
from pathlib import Path
""" Script to write fasta files with the tips left in the tree after pruning
Input   snakemake.input.seq_dict : pickled {tip label: Bio.SeqRecord} dictionary
        snakemake.input.tree     : pruned Newick tree
Output  snakemake.output.fasta   : sequences of the retained tips
"""
with open(snakemake.input.seq_dict, "rb") as f:
  prot_seq = pickle.load(f)
  
nw = snakemake.input.read_text().strip()
if not nw.endswith(";"):
    nw += ";"
t = PhyloTree(nw)
all_leaves = set(t.leaf_names())
records = [prot_seq[leaf] for leaf in all_leaves if leaf in prot_seq]
if records:
    output_path = Path(snakemake.output.fasta)
    with open(output_path, 'w') as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    print(f"{len(records)} sequences > {output_path}")
