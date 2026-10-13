"""
domain_mapping.py -- draw a subtree alongside the domain architecture of its tips.

Purpose
-------
Topology alone does not establish that the members of a candidate clade share
the architecture of the family under study: a clade can be recovered on the
strength of one conserved catalytic domain while its members differ in their
accessory domains. This script renders a gene (sub)tree next to a scaled
schematic of each tip protein, with annotated domains coloured by accession, so
that topology and architecture can be inspected together.

Input
-----
--tree    Newick file; tip labels are expected in the form <protein_id>@<taxon>,
          the protein identifier being the key used in the annotation table.
--tsv     InterProScan tab-separated output (protein id, _, sequence length,
          analysis/member database, signature accession, _, start, stop, ...).
--db      Member database to plot, e.g. Pfam (default) or SUPERFAMILY; rows from
          other analyses are read for sequence length only and not drawn.

Output
------
--outfile A raster image (PNG) containing the cladogram, tip labels, and one
          domain track per tip with a legend of the plotted accessions.

Procedure
---------
1. Parse the InterProScan table into {protein: {length, [(accession, start, stop)]}},
   taking the sequence length as the maximum reported across that protein's rows.
2. Compute a rectangular cladogram layout: tips are spaced evenly on the vertical
   axis in their Newick order and placed at x = 0; each internal node is placed at
   x = -(node depth) and at the mean vertical position of its children.
3. Assign one colour per domain accession from a fixed qualitative palette,
   cycling if the number of accessions exceeds the palette.
4. Draw the tree as elbow connectors (a vertical segment spanning the children's
   positions, then a horizontal segment to each child).
5. For each tip, write the label and draw a grey bar of fixed width representing
   the full-length protein, with each annotated domain overlaid at its
   proportional start and end coordinates.

Notes and caveats
-----------------
- The layout is a cladogram: horizontal positions encode topological depth, not
  branch length, so no evolutionary distance should be read off the figure.
- Domain tracks are normalised to a constant bar width, so protein lengths are
  not comparable between tips; only the relative position of domains within a
  protein is to scale.
- Tips absent from the annotation table are drawn without a track and without
  warning.
- Only the first tree in the file is read (the text is split at the first
  semicolon).
"""

import os
import sys
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from ete4 import Tree

def parse_args():
    """Define and parse the command-line interface."""
    parser = argparse.ArgumentParser(
        description="Mapping domains onto the subtrees"
    )
    parser.add_argument("--tree",       required=True, help="Path to newick tree file")
    parser.add_argument("--tsv",        required=True, help="Path to interpro tsv file")
    parser.add_argument("--outfile",    required=True, help="Path to output file (PNG)")
    parser.add_argument("--db",type=str, default='Pfam',   help="Domain databses to be considered (default: Pfam)")
    return parser.parse_args()

def parse_interpro_tsv(tsv_path, dom_db):
    """Index an InterProScan table by protein.

    Parameters
    ----------
    tsv_path : str
        Path to InterProScan tabular output. A missing file yields an empty index
        rather than an error, so the tree can still be drawn without annotation.
    dom_db : str
        Member-database name to retain (field 4 of the table).

    Returns
    -------
    defaultdict
        {protein_id: {"length": int, "domains": [(accession, start, stop), ...]}}.
        Rows from other member databases still contribute their sequence length,
        which keeps the length estimate available even when the requested database
        has no hit for that protein.
    """
    data = defaultdict(lambda: {'length': 0, 'domains': []})
    if not os.path.exists(tsv_path):
        return data
    with open(tsv_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            protein_id = parts[0]
            seq_len = int(parts[2])
            analysis = parts[3]
            pfam_id = parts[4]
            start = int(parts[6])
            stop = int(parts[7])
            
            # Sequence length is reported on every row; take the maximum defensively
            # in case of inconsistent records.
            data[protein_id]['length'] = max(data[protein_id]['length'], seq_len)
            if analysis == dom_db:
                data[protein_id]['domains'].append((pfam_id, start, stop))
        return data

def compute_rectangular_layout(tree):
    """Assign (x, y) plotting coordinates to every node of the tree.

    Tips are laid out at x = 0 and at consecutive integer y positions in Newick
    order. Internal nodes are placed at x = -depth and at the mean y of their
    children, which is computed by traversing the tree in reverse level order so
    that all children are positioned before their parent.

    Returns
    -------
    dict
        {node: (x, y)}.
    """
    coords = {}
    leaves = list(tree.leaves())
    
    for i, leaf in enumerate(leaves):
        coords[leaf] = (0, i)
    depths = {}
    # Topological depth (number of edges from the root), used as the x coordinate.
    def calc_depth(node, d=0):
        depths[node] = d
        for child in node.children:
            calc_depth(child, d + 1)
    
    calc_depth(tree)
    for node in reversed(list(tree.traverse())):
        if node in coords:
            continue
        
        children_coords = [coords[c] for c in node.children]
        avg_y = sum(c[1] for c in children_coords) / len(children_coords)
        coords[node] = (-depths[node], avg_y)
        
    return coords

def main():
    """Render the annotated tree figure."""
    args = parse_args()
    with open(args.tree, 'r') as f:
        tree_str = f.read().split(';')[0].strip() + ';'
    t = Tree(tree_str)
    domains_data = parse_interpro_tsv(args.tsv, args.db)
    coords = compute_rectangular_layout(t)
    leaves = list(t.leaves())
    fig, ax = plt.subplots(figsize=(12, len(leaves) * 0.3))
    ax.set_axis_off()
    # Qualitative palette for domain accessions; the figure height scales with the
    # number of tips so that label spacing stays constant.
    palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999', '#000000']
    all_pfams = set()
    for prot_id in domains_data:
        for pfam, _, _ in domains_data[prot_id]['domains']:
            all_pfams.add(pfam)
    all_pfams = sorted(list(all_pfams))
    pfam_color_map = {pfam: palette[i % len(palette)] for i, pfam in enumerate(all_pfams)}
    # Step 4: elbow connectors. For each internal node, one vertical segment through
    # the span of its children, then one horizontal segment out to each child.
    for node in t.traverse():
        if not node.children:
            continue
        
        x1, y1 = coords[node]
        for child in node.children:
            x2, y2 = coords[child]
            ax.plot([x1, x1], [y1, y2], color='black', linewidth=1, zorder=1)
            ax.plot([x1, x2], [y2, y2], color='black', linewidth=1, zorder=1)
    # Step 5: tip labels and domain tracks. The protein identifier is recovered from
    # the part of the tip label preceding '@', which is the key used by InterProScan.
    for leaf in leaves:
        x, y = coords[leaf]
        label = leaf.name
        prot_id = label.split('@')[0]
        
        ax.text(x + 0.1, y, label, fontsize=8, va='center', ha='left')
        if prot_id in domains_data:
            info = domains_data[prot_id]
            prot_len = info['length']
            
            if prot_len > 0:
                # Fixed horizontal offset and width for every protein bar: domain
                # coordinates are rescaled onto this constant-width track, so tracks
                # are internally to scale but not comparable between proteins.
                bar_start = x + 5 
                bar_width = 3
                ax.add_patch(patches.Rectangle((bar_start, y - 0.1), bar_width, 0.2, color='lightgrey', zorder=1))
                for pfam, start, stop in info['domains']:
                    rel_start = bar_start + (start / prot_len) * bar_width
                    rel_stop = bar_start + (stop / prot_len) * bar_width
                    color = pfam_color_map.get(pfam, '#cccccc')
                    ax.add_patch(patches.Rectangle((rel_start, y - 0.1), rel_stop - rel_start, 0.2, color=color, zorder=2))
    
    if all_pfams:
        legend_patches = [patches.Patch(color=pfam_color_map[pfam], label=pfam) for pfam in all_pfams]
        ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., title=f"{args.db} Domains")
    
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=96, bbox_inches='tight')

if __name__ == "__main__":
    main()
