"""
colour_subtree.order.py to annotate and colour a selected subtree for figure
preparation.

Input
-----
--tree        Newick subtree (one tree, tip labels as database accessions).
--lineage     NCBI fullnamelineage.dmp (taxid | name | full lineage |).
--asg_lineage Curated mapping of Asgard/Thermoproteota record identifiers to
              taxon name and phylum (tab-separated).
--taxa_map    Table used to attach a GTDB-style label to focal tips; must
              contain the columns 'extracted_classification' and 'short_taxon'.
--rep_db      SQLite database (table `dictionary`, columns key/value) for the ncbi prot2acessiontaxID file, mapping
              cluster representative accessions to taxon names.

Output
------
--output      NEXUS tree with per-tip colour attributes and FigTree settings.
<output>.taxa.tsv  Per-tip table of taxa domain

"""

import argparse
import sqlite3
import re
import random
import csv
import sys
from pathlib import Path
# The shared helpers (get_species_name, to_nexus) live two directories above this
# script; prepend that location to the module search path before importing them.
target_dir = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(target_dir))

from focal_subtree_build import get_species_name, to_nexus
from ete4 import Tree


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
    parser = argparse.ArgumentParser(description="Generate colored/annotated tree + tip metadata for subtree")
    parser.add_argument("--tree", required=True, help="Path to Newick subtree file")
    parser.add_argument("--lineage", default="/local/two/vineeth/DNArep_asCOG-nr/gene_trees/focal_subTree/fullnamelineage.dmp", 
                        help="Path to lineage file (| delimited)")
    parser.add_argument("--asg_lineage", default="/local/two/vineeth/DNArep_asCOG-nr/gene_trees/focal_subTree/cleaned_asgard+thermo_mapping.tsv", 
                        help="Path to asgard mapping file")
    parser.add_argument("--taxa_map", default="/local/two/vineeth/DNArep_asCOG-nr/gene_trees/focal_subTree/extracted_taxa.tsv",
                         help="Path to taxa map file, used to label Promethearchaeati tips")
    parser.add_argument("--output", required=True, help="Path to write the tree")
    parser.add_argument("--rep_db", default="/local/two/vineeth/DNArep_asCOG-nr/gene_trees/subTrees/rep_dict.db",
                         help="Path to the representative-name SQLite lookup db")
    return parser.parse_args()

def build_lineage_dicts(lineage_file, map_file):
    """Build the taxonomy lookup dictionaries.
    """
    lineage_dict = {}
    name_to_taxid = {}
    asg_map_dict = {}
    with open(lineage_file, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 3:
                tax_id = parts[0]
                taxon_name = parts[1].replace(" ", "_")
                lineage = parts[2]
                lineage_dict[tax_id] = {
                    'name': taxon_name,
                    'lineage': lineage,
                    'domain': get_domain(lineage),
                    'phylum': get_phylum(lineage),
                    'order': get_order(lineage)
                }
                name_to_taxid[taxon_name] = tax_id
    with open(map_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                record_id = parts[0]
                taxon_name = parts[1]
                asg_map_dict[record_id] = {
                    'name': taxon_name.replace(" ", "_"),
                    'phylum': parts[2]
                }
    return lineage_dict, name_to_taxid, asg_map_dict

def get_domain(lineage):
    """Return the domain of life by substring search of the lineage string."""
    for domain in ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']:
        if domain in lineage:
            return domain
    return 'Unknown'

def get_phylum(lineage):
    """Return the third semicolon-delimited lineage field as the phylum based on the ncbi accession.
    """
    parts = lineage.split(';')
    return parts[2].strip() if len(parts) >= 3 else 'Unknown'
    
def get_order(lineage):
    parts = lineage.split(';')
    return parts[3].strip() if len(parts) >= 4 else 'Unknown'

def get_tax_id(name, tax_id_dict):
    """Resolve the taxon name embedded in an <accession>@<taxon> label to a taxid."""
    parts = name.split('@')
    if len(parts) < 2:
        print(f"Warning: skipping {name}")
        return None
    taxon_name = parts[1]
    return tax_id_dict.get(taxon_name)
    
gtdb_pattern = re.compile(r'_sp\d{6,}$')
def build_extracted_taxa_dict(tsv_path):
    """Collapse the taxa map to one short label per extracted classification.
    """
    from collections import Counter, defaultdict
    counts = defaultdict(Counter)
    with open(tsv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'extracted_classification', 'short_taxon'}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"Expected columns {required} in {tsv_path}, got: {reader.fieldnames}")
        for row in reader:
            key = (row.get('extracted_classification') or '').strip()
            short_taxon = (row.get('short_taxon') or '').strip()
            if key and short_taxon:
                counts[key][short_taxon] += 1
    return {key: counter.most_common(1)[0][0] for key, counter in counts.items()}

def match_extracted_taxon(taxon_name, taxa_dict):
    """Look up a focal taxon in the taxa map, finding the GTDB suffixes.
    """
    query = taxon_name.replace('_', ' ').strip()
    if query in taxa_dict:
        return taxa_dict[query]
    stripped = gtdb_pattern.sub('', taxon_name).replace('_', ' ').strip()
    if stripped and stripped in taxa_dict:
        return taxa_dict[stripped]
    return None

def color_and_export_subtree(tree, lineage_dict, name_to_taxid, asg_map_dict, accession_db,
                              output_path, taxa_dict=None, order_to_color=None, color_state=None):
    """Rewrite, classify and colour every tip, then export tree and metadata.
    """
    if order_to_color is None:
        order_to_color = {
            'Promethearchaeati': '#C3B1E1',
            'Unknown': '#708090'
        }
    if color_state is None:
        color_state = {'index': 0}
 
    color_palette = [
        '#143349', '#DAA464', '#ab665a', '#049558', '#ffbb8b',
        '#56454b', '#FF5FCF', '#1f184e', '#F25912', '#99566e',
        '#827717', '#607d8b', '#9a6324', '#566d99', '#215E61',
        '#aaffc3', '#808000', '#5C3E94', '#767F9E', '#808080'
    ]
 
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    name_to_phylum = {data['name']: data['phylum'] for data in asg_map_dict.values()}
    metadata_path = output_path.with_suffix('.taxa.tsv')
    with open(metadata_path, 'w') as out_file:
        out_file.write("Taxon\tPhylum\tDomain\tGTDB_Taxon\n")
        
        for leaf in tree.leaves():
            if not leaf.name:
                leaf.detach()
                continue
            leaf.name = get_species_name(leaf.name, lineage_dict, accession_db, asg_map_dict)
            leaf_name = leaf.name
            if leaf_name.startswith('ref|'):
                taxon_name = leaf_name
                domain = 'Reference'
                phylum = 'Reference'
                order = 'Reference'
                matched_taxon = ''
                extra_tag = '{Reference}'
                leaf.name = leaf.name + extra_tag
                color_key = 'Reference'
                if color_key not in order_to_color:
                    if color_state['index'] < len(color_palette):
                        order_to_color[color_key] = color_palette[color_state['index']]
                        color_state['index'] += 1
                    else:
                        order_to_color[color_key] = "#{:06x}".format(random.randint(0, 0xFFFFFF))
                color = order_to_color[color_key]
                leaf.add_prop('color', color)
                out_file.write(f"{taxon_name}\t{order}\t{phylum}\t{domain}\t{matched_taxon}\n")
                continue

            if '@' not in leaf_name:
                print(f"Warning: skipping {leaf_name}")
                taxon_name = leaf_name
                taxon_id = None
            else:
                taxon_name = leaf_name.split('@')[1]
                taxon_id = get_tax_id(leaf_name, name_to_taxid)

            domain = 'Unknown'
            phylum = 'Unknown'
            order = 'Unknown'
            found_taxid = None
            if taxon_name in name_to_phylum:
                domain = 'Archaea'
                exp_phylum = name_to_phylum[taxon_name]
                phylum = 'Promethearchaeati' if exp_phylum == 'Asgardarchaeota' else exp_phylum
            elif taxon_id in lineage_dict:
                found_taxid = taxon_id
            if found_taxid and found_taxid in lineage_dict:
                domain = lineage_dict[found_taxid].get('domain', 'Unknown')
                phylum = lineage_dict[found_taxid].get('phylum', 'Unknown')
                order = lineage_dict[found_taxid].get('order', 'Unknown')
            if not domain:
                domain = 'Unknown'
            if not phylum or "incertae sedis" in phylum.lower() or "unclassified" in phylum.lower():
                phylum = "Unknown"
            if phylum == 'Promethearchaeati':
                matched_taxon = ''
                if taxa_dict:
                    matched_taxon = match_extracted_taxon(taxon_name, taxa_dict) or ''
                    if not matched_taxon:
                        print(f"Warning: no taxa match found for Promethearchaeati tip '{taxon_name}'")
                extra_tag = '{Phylum: ' + phylum
                if matched_taxon:
                    safe_matched_taxon = matched_taxon.replace(';', '|')
                    extra_tag += '//' + safe_matched_taxon
                extra_tag += '}'
            else:
                extra_tag = f"{{Phylum: {phylum}//{order}}}"
                matched_taxon = '' 
            leaf.name = leaf.name + extra_tag
            color_key = phylum
            if color_key not in order_to_color:
                if color_state['index'] < len(color_palette):
                    order_to_color[color_key] = color_palette[color_state['index']]
                    color_state['index'] += 1
                else:
                    order_to_color[color_key] = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            color = order_to_color[color_key]
            leaf.add_prop('color', color)
            out_file.write(f"{taxon_name}\t{order}\t{phylum}\t{domain}\t{matched_taxon}\n")
    to_nexus(tree, output_path)
    return order_to_color


def main():
    """Load the subtree and reference data, annotate, and write the outputs."""
    args = parse_args()
    accession_db = LargeDictionary(args.rep_db)
    
    tree_str = Path(args.tree).read_text().strip()
    if not tree_str.endswith(";"):
        tree_str += ";"
    t = Tree(tree_str)

    lineage_dict, name_to_taxid, asg_map_dict = build_lineage_dicts(args.lineage, args.asg_lineage)
    taxa_dict = build_extracted_taxa_dict(args.taxa_map)
 
    color_and_export_subtree(
        t, lineage_dict, name_to_taxid, asg_map_dict, accession_db,
        Path(args.output), taxa_dict=taxa_dict
    )
    accession_db.close()


if __name__ == "__main__":
    main()