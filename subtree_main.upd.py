"""
subtree_main.upd.py to identify focal-clade (Asgard) subtree selection from a
single gene tree.

Input
-----
--tree        Newick gene tree with database accessions as tip labels.
--lineage     NCBI fullnamelineage.dmp (taxid | name | full lineage |).
--asg_lineage Curated Asgard record-to-taxon mapping (tab-separated).
--rep_db      SQLite accession-to-taxon lookup (table `dictionary`).
--out         Output directory; its basename is used as the file prefix.

Output (written to --out)
-------------------------
<prefix>_subtree_<rank>.nwk        selected subtree, Newick
<prefix>_subtree_<rank>_tip_ids.txt  accessions of its tips
<prefix>_subtree_<rank>.nex        the same subtree, coloured for FigTree
<prefix>_subtree_<rank>_taxa.txt   per-tip taxonomy of that subtree
<prefix>_report.json               run parameters and every candidate clade
<prefix>_candidate_overview.tsv    one row per candidate clade for manual review
<prefix>_subtree_<rank>.fasta      tip sequences, only if MASTER_FASTA is set

"""

import argparse
from os.path import basename
from pathlib import Path
import sqlite3
from subtree_selection import (
    build_lineage_dicts,
    build_focal_list,
    find_all_focal_taxa_in_tree,
    find_monophyletic_focal_clades,
    get_trusted_family_region,
    select_monophyletic_asgard_subtrees,
    write_full_report,
    write_candidate_overview_tsv,
    color_and_export_subtree,
    write_tip_ids,
    extract_fasta_by_ids,
)
# ---------------------------------------------------------------------------
# Run parameters.
# ---------------------------------------------------------------------------
TARGET_FAMILY = None          
MASTER_FASTA = None             
MIN_MONOPHYLETIC_SIZE = 4       

MAX_SISTER_TAXA = 300           
MAX_SUBTREE_SIZE = 450          
MIN_TRUSTED_REFS = 1            
COMPOSITION_THRESHOLDS = (75, 90, 100)
SISTER_ABSORB_THRESHOLD = 90.0
FOCAL_PHYLUM = 'Promethearchaeati'


class LargeDictionary:
    def __init__(self, db_path="rep_dict.db"):
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()

    def lookup(self, key):
        self.cursor.execute('SELECT value FROM dictionary WHERE key = ?', (key,))
        result = self.cursor.fetchone()
        return result[0] if result else None

    def close(self):
        self.conn.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Asgard subtree selection")
    parser.add_argument("--tree", required=True, help="Path to newick tree file")
    parser.add_argument("--lineage", required=True, help="Path to lineage file (| delimited)")
    parser.add_argument("--asg_lineage", required=True, help="Path to asgard mapping file")
    parser.add_argument("--out", required=True, help="Directory for output files")
    parser.add_argument("--rep_db", required=True, help="Path to the representative-name SQLite lookup db")
    return parser.parse_args()


def main():
    """Execute the eight steps of the selection procedure for one gene tree."""
    args = parse_args()
    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)
    ncbi_dict = LargeDictionary(args.rep_db)
    prefix_base = basename(str(output_dir))

    # Step 1: taxonomy lookups, then the focal taxon list
    lineage_dict, name_to_taxid, asg_dict = build_lineage_dicts(args.lineage, args.asg_lineage)
    focal_list = build_focal_list(lineage_dict, asg_dict)
    if not focal_list:
        return

    # Step 2: relabel tips as <accession>@<taxon> and identify the focal tips.
    t, focal_taxa_in_tree = find_all_focal_taxa_in_tree(
        args.tree, focal_list, lineage_dict, ncbi_dict, asg_dict
    )
    if not focal_taxa_in_tree:
        print(f'No asgards found in tree: {args.tree}')
        return

    # Step 3: optional restriction to a trusted family region.
    trusted_leaf_names = None
    trusted_clusters_report = []
    if TARGET_FAMILY:
        trusted_leaf_names, trusted_clusters_report = get_trusted_family_region(
            t, TARGET_FAMILY, min_refs=MIN_TRUSTED_REFS
        )
        if trusted_leaf_names is None:
            print(f"WARNING: no reference cluster for '{TARGET_FAMILY}' had >= "
                  f"{MIN_TRUSTED_REFS} contamination-free members -- proceeding "
                  f"WITHOUT a trusted-region restriction. See 'trusted_family_clusters' "
                  f"in the report and consider lowering MIN_TRUSTED_REFS or curating "
                  f"references manually.")
        else:
            print(f"Trusted '{TARGET_FAMILY}' region: {len(trusted_leaf_names)} leaves.")

    # Step 4: maximal focal-only clades above the size threshold, largest first.
    monophyletic_clades = find_monophyletic_focal_clades(
        t, focal_taxa_in_tree, min_size=MIN_MONOPHYLETIC_SIZE,
        trusted_leaf_names=trusted_leaf_names
    )
    if not monophyletic_clades:
        print(f'No monophyletic asgard clades with >= {MIN_MONOPHYLETIC_SIZE} members '
              f'found in tree: {args.tree}'
              + (' within the trusted family region' if trusted_leaf_names is not None else ''))
        return

    focal_node_dict = {i: node for i, node in enumerate(monophyletic_clades)}

    # Step 5: context extraction, downsampling, composition summary and transfer-type
    selected, all_candidates = select_monophyletic_asgard_subtrees(
        t, focal_node_dict, focal_list, lineage_dict, name_to_taxid, asg_dict,
        max_subtree_size=MAX_SUBTREE_SIZE,
        max_sister_taxa=MAX_SISTER_TAXA,
        composition_thresholds=COMPOSITION_THRESHOLDS,
        sister_absorb_threshold=SISTER_ABSORB_THRESHOLD,
        focal_phylum=FOCAL_PHYLUM,
    )
    phylum_to_color = {}
    color_state = {'index': 0}
    id_to_fasta_paths = {}

    # Step 6: write Newick, tip identifiers and coloured NEXUS for each subtree.
    for entry in selected:
        prefix = f"{prefix_base}_subtree_{entry['subtree_rank']}"
        subtree_file = output_dir / f"{prefix}.nwk"
        entry['subtree'].write(outfile=str(subtree_file))

        ids_file = output_dir / f"{prefix}_tip_ids.txt"
        tip_ids = write_tip_ids(entry['subtree'], ids_file)
        if MASTER_FASTA:
            fasta_out = str(output_dir / f"{prefix}.fasta")
            for tip_id in tip_ids:
                id_to_fasta_paths.setdefault(tip_id, []).append(fasta_out)
        phylum_to_color = color_and_export_subtree(
            entry['subtree'], lineage_dict, name_to_taxid, asg_dict, ncbi_dict,
            output_dir, prefix, phylum_to_color=phylum_to_color, color_state=color_state
        )

    # Step 7: single pass over the master FASTA.
    if MASTER_FASTA and id_to_fasta_paths:
        found, missing = extract_fasta_by_ids(MASTER_FASTA, id_to_fasta_paths)
        print(f"FASTA extraction: {len(found)}/{len(id_to_fasta_paths)} tip IDs matched in {MASTER_FASTA}")
        if missing:
            print(f"WARNING: {len(missing)} tip IDs not found in master FASTA (showing up to 10): "
                  f"{list(missing)[:10]}")

    # Step 8: the JSON report records the parameters alongside every candidate
    # the TSV report gives one row per candidate for manual curation.
    run_info = {
        "tree": args.tree,
        "target_family": TARGET_FAMILY,
        "min_trusted_refs": MIN_TRUSTED_REFS,
        "min_monophyletic_size": MIN_MONOPHYLETIC_SIZE,
        "max_sister_taxa": MAX_SISTER_TAXA,
        "max_subtree_size": MAX_SUBTREE_SIZE,
        "composition_thresholds": COMPOSITION_THRESHOLDS,
        "sister_absorb_threshold": SISTER_ABSORB_THRESHOLD,
        "n_focal_taxa_in_tree": len(focal_taxa_in_tree),
        "n_candidates_found": len(all_candidates),
        "n_subtrees_selected": len(selected),
    }
    report_path = output_dir / f"{prefix_base}_report.json"
    write_full_report(report_path, run_info, trusted_clusters_report, all_candidates)
    overview_path = output_dir / f"{prefix_base}_candidate_overview.tsv"
    write_candidate_overview_tsv(
        all_candidates, overview_path,
        lineage_dict, name_to_taxid, asg_dict,
        family_type=TARGET_FAMILY if TARGET_FAMILY else 'NA',
    )
    ncbi_dict.close()


if __name__ == "__main__":
    main()