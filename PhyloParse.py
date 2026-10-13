"""
Prunes out-paralogues and in-paralogues from a eukaryotic gene tree, so that each
tree approximates one LECA orthologous group.

Input   snakemake.input.tree  : Newick gene tree, tips as <prefix>_<species>_<seq id>
        snakemake.input.guide : species guide, columns 'Name' and 'Group'
Output  snakemake.output.tree : pruned Newick tree

"""
from ete4 import PhyloTree
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import re

input_file = snakemake.input.tree
output_file = snakemake.output.tree

def load_supergroups(guide_path):
    df = pd.read_csv(guide_path, sep='\t')
    species_names = df['Name'].str.split('_').str[1:-1].str.join('_')
    return dict(zip(species_names, df['Group']))
SUPERGROUP_MAP = load_supergroups(snakemake.input.guide)
def species_name(name):
    parts = name.split('_')
    return "_".join(parts[1:-1])
def get_supergroup(species):
    return SUPERGROUP_MAP.get(species)
def leca_node(seqs):
    groups = {get_supergroup(species_name(s)) for s in seqs if species_name(s)}
    return 'Opimoda' in groups and 'Diphoda' in groups
def ancestral_dup(node,dup_nodes):
    return not any(anc in dup_nodes for anc in node.ancestors())
def tree_prune(input, output_tree):
    input_path = Path(input)
    output_tree_path = Path(output_tree)
    output_tree_path.parent.mkdir(parents=True, exist_ok=True)
    nw = input_path.read_text().strip()
    if not nw.endswith(";"):
        nw += ";"
    nw = re.sub(r'\)([0-9.]+):', r'):', nw)
    t = PhyloTree(nw)
    t.set_species_naming_function(species_name)
    num_leaves = len(list(t.leaves()))
    midpoint = t.get_midpoint_outgroup()
    if midpoint != t:
        t.set_outgroup(midpoint)
    for node in t.traverse():
        if len(node.children) > 2:
            node.resolve_polytomy()
    all_events = list(t.get_descendant_evol_events())
    dup_nodes = {ev.node for ev in all_events if ev.etype == 'D'}
    out_paralog = set()
    for ev in all_events:
        if ev.etype == 'D':
            if leca_node(ev.in_seqs) and leca_node(ev.out_seqs):
                if ancestral_dup(ev.node, dup_nodes):
                    dropped_clade = set(min(ev.in_seqs, ev.out_seqs, key=len))
                    branches_removed = out_paralog | dropped_clade
                    if len(branches_removed) >= (num_leaves*0.75):
                        continue
                    out_paralog.update(dropped_clade)
    if out_paralog:
        leaves_to_keep = [l.name for l in t.leaves() if l.name not in out_paralog]
        if leaves_to_keep:
            t.prune(leaves_to_keep)
        else:
            print(f"Skipping {input_path.name}: No leaves left after out-paralog pruning.")
    pruning_active = True
    while pruning_active:
        pruning_active = False
        all_events = list(t.get_descendant_evol_events())
        for ev in all_events:
            if ev.etype == 'D':
                in_species = {species_name(s) for s in ev.in_seqs}
                out_species = {species_name(s) for s in ev.out_seqs}
                overlap = in_species & out_species
                if overlap:
                    smaller_clade = min(in_species, out_species, key=len)
                    overlap_ratio = len(overlap) / len(smaller_clade)
                    if overlap_ratio >= 0.5:
                        target = ev.in_seqs if len(in_species) <= len(out_species) else ev.out_seqs
                        leaves_to_prune = {s for s in target if species_name(s) in overlap}
                        if leaves_to_prune:
                            remaining_leaves = [l.name for l in t.leaves() if l.name not in leaves_to_prune]
                            if remaining_leaves:
                                t.prune(remaining_leaves)
                                pruning_active = True
                                break
        if len(list(t.leaves())) <= 1:
            break
    t.write(outfile=str(output_tree_path))
    print(f"Processed {input_path.name} > {output_tree_path.name}")

tree_prune(input_file, output_file)
