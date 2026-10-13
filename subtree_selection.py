"""
subtree_selection.py, the backend script for subtree_main.upd, colour_subtree.order for the delimitation, extraction and classification of focal
(Asgard) clades within a gene tree.
"""

from ete4 import Tree
from pathlib import Path
from collections import Counter
import re
import json
def build_lineage_dicts(lineage_file, map_file):
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
                    'phylum': get_phylum(lineage)
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
    for domain in ['Bacteria', 'Archaea', 'Eukaryota', 'Viruses']:
        if domain in lineage:
            return domain
    return 'Unknown'
def get_phylum(lineage):

    parts = lineage.split(';')
    return parts[2].strip() if len(parts) >= 3 else 'Unknown'
def get_tax_id(name, tax_id_dict):
    parts = name.split('@')
    try:
        taxon_name = parts[1]
        tax_id = tax_id_dict.get(taxon_name)
        return tax_id
    except IndexError:
        print(f"IndexError triggered by label -> {name}")
        raise
def build_focal_list(lineage_dict, asg_map_dict):
    """
    Returns a list of (identifier, taxon name) pairs comprising the curated records
    assigned to Asgardarchaeota and every NCBI taxon whose lineage contains
    Promethearchaeati. Identifiers and names are matched separately downstream, so
    a tip resolvable by either route is recognised as focal.
    """
    focal_list = []
    for record_id, data in asg_map_dict.items():
        if data['phylum'] == 'Asgardarchaeota':
            focal_list.append((record_id, data['name']))
    for tax_id, info in lineage_dict.items():
        if 'Promethearchaeati' in info['lineage']:
            focal_list.append((tax_id, info['name']))
    return focal_list
def get_species_name(label, lineage_dict, rep_dict, asg_map_dict):
    """Append the source organism to a tip label, returning <label>@<taxon>.
    """
    if label.startswith('ref|'):
        parts = label.split('|')
        species = parts[2] if len(parts) >= 3 else label
        return f"{label}@{species}"
    asg_pattern = re.compile(r'^(.*?)_+\d+$')
    taxid = label.split('.')[0]
    if taxid in lineage_dict:
        new_label = label + '@' + lineage_dict[taxid]['name']
        return new_label
    taxon_name = rep_dict.lookup(label)
    if taxon_name:
        taxon_name = taxon_name.replace(" ", "_")
        new_label = f"{label}@{taxon_name}"
        return new_label
    asg_match = asg_pattern.search(label)
    if asg_match:
        asg_id = asg_match.group(1)
        asgard_taxon = asg_map_dict.get(asg_id)
        if asgard_taxon:
            new_label = f"{label}@{asgard_taxon['name']}"
            return new_label
    return label
def find_all_focal_taxa_in_tree(tree_file, focal_list, lineage_dict, rep_dict, asg_map_dict):
    """Read a gene tree, relabel its tips, and list the focal tips.
    """
    tree_path = Path(tree_file)
    nw = tree_path.read_text().strip()
    if not nw.endswith(";"):
        nw += ";"
    t = Tree(nw)
    for leaf in t.leaves():
        leaf.name = get_species_name(leaf.name, lineage_dict, rep_dict, asg_map_dict)
    focal_taxa_in_tree = []
    focal_ids = {tax_id for tax_id, name in focal_list}
    focal_names = {name for tax_id, name in focal_list}
    for node in t.traverse():
        if node.is_leaf:
            parts = node.name.split('@')
            if len(parts) >= 2:
                taxon_name = parts[1]
                id_part = parts[0]
                ncbi_id = id_part.split('.')[0]
                if taxon_name in focal_names or ncbi_id in focal_ids:
                    focal_taxa_in_tree.append(node.name)
    return t, focal_taxa_in_tree


def find_monophyletic_focal_clades(t, focal_taxa_in_tree, min_size=4, trusted_leaf_names=None):
    """Delimit the clades composed exclusively of focal tips.
    Returns the qualifying nodes, largest first.
    """
    focal_set = set(focal_taxa_in_tree)
    focal_only_nodes = []
    for node in t.traverse():
        leaf_names = set(leaf.name for leaf in node.leaves())
        if leaf_names and leaf_names.issubset(focal_set):
            focal_only_nodes.append(node)

    focal_only_ids = set(id(n) for n in focal_only_nodes)
    unique_nodes = [
        n for n in focal_only_nodes
        if not any(id(anc) in focal_only_ids for anc in n.ancestors())
    ]

    qualifying = [n for n in unique_nodes if len(list(n.leaves())) >= min_size]

    if trusted_leaf_names is not None:
        qualifying = [
            n for n in qualifying
            if set(l.name for l in n.leaves()).issubset(trusted_leaf_names)
        ]

    qualifying.sort(key=lambda n: len(list(n.leaves())), reverse=True)
    return qualifying


def find_trusted_family_clusters(t, target_family, ref_prefix='ref|'):
    """Group the curated references of one family into contamination-free clusters.
    """
    all_refs = [l for l in t.leaves() if l.name.startswith(ref_prefix)]
    target_refs = [r for r in all_refs if r.name.split('|')[1] == target_family]
    other_names = set(r.name for r in all_refs if r.name.split('|')[1] != target_family)

    ref_by_id = {id(r): r for r in target_refs}
    unassigned = set(ref_by_id)
    clusters = []
    while unassigned:
        seed_id = next(iter(unassigned))
        cluster = [ref_by_id[seed_id]]
        unassigned.discard(seed_id)
        progress = True
        while progress:
            progress = False
            for cand_id in list(unassigned):
                trial = cluster + [ref_by_id[cand_id]]
                mrca = t.common_ancestor(trial) if len(trial) > 1 else trial[0]
                leafnames = set(l.name for l in mrca.leaves())
                if not (leafnames & other_names):
                    cluster.append(ref_by_id[cand_id])
                    unassigned.discard(cand_id)
                    progress = True
        clusters.append(cluster)

    def _clade(cluster):
        return t.common_ancestor(cluster) if len(cluster) > 1 else cluster[0]

    clusters.sort(key=lambda c: len(list(_clade(c).leaves())), reverse=True)
    return [(_clade(c), c) for c in clusters]


def get_trusted_family_region(t, target_family, min_refs=2):
    """Select the trusted region of the tree for one family.
    """
    clusters = find_trusted_family_clusters(t, target_family)

    clusters_report = []
    trusted_leaf_names = None
    chosen = False
    for i, (clade, refs) in enumerate(clusters):
        is_trusted = (not chosen) and len(refs) >= min_refs
        if is_trusted:
            chosen = True
            trusted_leaf_names = set(l.name for l in clade.leaves())
        clusters_report.append({
            "cluster_rank": i + 1,
            "n_refs": len(refs),
            "clade_size": len(list(clade.leaves())),
            "trusted": is_trusted,
            "ref_names": [r.name for r in refs],
        })

    return trusted_leaf_names, clusters_report


def get_subtree_with_context(focal_node):
    """Return the subtree used for origin inference report, with its context nodes.
    """
    parent = focal_node.up
    if parent is None:
        return focal_node, [], []
    sister_nodes = [c for c in parent.children if c is not focal_node]
    grandparent = parent.up
    if grandparent is None:
        return parent, sister_nodes, []
    outgroup_nodes = [c for c in grandparent.children if c is not parent]
    return grandparent, sister_nodes, outgroup_nodes


def _clade_focal_fraction(node, lineage_dict, name_to_taxid, name_to_phylum, focal_phylum):
    """Percentage of a clade's tips assigned to the focal phylum."""
    leaves = list(node.leaves())
    if not leaves:
        return 0.0
    n_match = 0
    for leaf in leaves:
        _, phylum = get_leaf_domain_phylum(leaf.name, lineage_dict, name_to_taxid, name_to_phylum)
        if phylum == focal_phylum:
            n_match += 1
    return 100.0 * n_match / len(leaves)


def expand_focal_clade(focal_node, lineage_dict, name_to_taxid, asg_map_dict,
                        threshold=90.0, focal_phylum='Promethearchaeati'):
    """Climbs from focal_node one node at a time, absorbing the parent only while all
    of the current node's siblings consist of at least `threshold` per cent focal
    tips. This recovers a focal clade that is interrupted by a small number of
    misplaced non-focal sequences.
    """
    name_to_phylum = {data['name']: data['phylum'] for data in asg_map_dict.values()}
    current = focal_node
    while True:
        parent = current.up
        if parent is None:
            break
        siblings = [c for c in parent.children if c is not current]
        if not siblings:
            break
        fractions = [
            _clade_focal_fraction(s, lineage_dict, name_to_taxid, name_to_phylum, focal_phylum)
            for s in siblings
        ]
        if all(f >= threshold for f in fractions):
            current = parent
        else:
            break
    return current


def absorb_dominant_sister_clades(focal_node, lineage_dict, name_to_taxid, asg_map_dict,
                                   target_phylum='Promethearchaeati', min_dominant_pct=0.0,
                                   max_iterations=50):
    """Absorb successive sister clades whose dominant phylum is the focal phylum.
    """
    name_to_phylum = {d['name']: d['phylum'] for d in asg_map_dict.values()}
    current = focal_node
    n_absorbed = 0
    absorbed_taxa = []

    for _ in range(max_iterations):
        parent = current.up
        if parent is None:
            break
        sister_nodes = [c for c in parent.children if c is not current]
        if not sister_nodes:
            break
        sister_leaf_names = [l.name for s in sister_nodes for l in s.leaves()]
        if not sister_leaf_names:
            break

        phyla = [
            get_leaf_domain_phylum(name, lineage_dict, name_to_taxid, name_to_phylum)[1]
            for name in sister_leaf_names
        ]
        dom_phylum, dom_pct = _dominant_value(phyla)

        if dom_phylum == target_phylum and dom_pct >= min_dominant_pct:
            absorbed_taxa.extend(sister_leaf_names)
            n_absorbed += 1
            current = parent
        else:
            break

    return current, n_absorbed, absorbed_taxa


def _find_N_bounded(boundary, leaf):
    """List the distances from one tip to its nearest relatives within a boundary.

    Implements the neighbour enumeration of the Treemmer downsampling heuristic
    (Menardo et al. 2018), restricted to the clade `boundary` so that pruning never
    considers tips outside the subtree being reduced. Sister tips are listed
    directly as a dictionary.
    """
    dlist = {}
    parent = leaf.up
    if parent is None or leaf is boundary:
        return dlist
    sister_count = 0
    niece_list = {}
    for child in parent.children:
        if child == leaf:
            continue
        if child.is_leaf:
            dlist[(leaf.name, child.name)] = boundary.get_distance(leaf, child)
            sister_count += 1
        else:
            for grandchild in child.children:
                if grandchild.is_leaf:
                    niece_list[(leaf.name, grandchild.name)] = boundary.get_distance(leaf, grandchild)
    if sister_count == 0 and len(niece_list) == 1:
        dlist.update(niece_list)
    if sister_count == 0 and parent is not boundary:
        grandparent = parent.up
        if grandparent is not None:
            aunt_dlist = {}
            for aunt in grandparent.children:
                if aunt == parent:
                    continue
                if aunt.is_leaf:
                    aunt_dlist[(leaf.name, aunt.name)] = boundary.get_distance(leaf, aunt)
            if len(aunt_dlist) == 1:
                dlist.update(aunt_dlist)
    return dlist


def treemmer_prune_to_n(boundary_node, target_n, protected_names=None):
    """Downsample a clade to 'n' tips by iterative nearest-neighbour pruning.
    At each iteration the closest pair of neighbouring tips is identified and one of
    the two removed, so that redundancy is lost before diversity: the tips discarded
    first are those whose phylogenetic position is most nearly duplicated by another
    tip. Within a pair, the tip with the shorter terminal branch is removed. The protected 
    names (focal tips, tips claimed by another subtree) are excluded from removal
    """
    protected_names = protected_names or set()
    while True:
        leaves = list(boundary_node.leaves())
        if len(leaves) <= target_n:
            break
        dlist = {}
        for leaf in leaves:
            dlist.update(_find_N_bounded(boundary_node, leaf))
        if not dlist:
            break
        # Prefer pairs that aren't both protected
        unprotected_pairs = {
            pair: dist for pair, dist in dlist.items()
            if not (pair[0] in protected_names and pair[1] in protected_names)
        }
        search_space = unprotected_pairs if unprotected_pairs else dlist
        pair = min(search_space, key=search_space.get)
        n1, n2 = pair
        node1, node2 = boundary_node[n1], boundary_node[n2]
        if n1 in protected_names and n2 not in protected_names:
            to_prune_name = n2
        elif n2 in protected_names and n1 not in protected_names:
            to_prune_name = n1
        else:
            to_prune_name = n1 if node1.dist < node2.dist else n2
        leaf_to_prune = boundary_node[to_prune_name]
        parent = leaf_to_prune.up
        leaf_to_prune.detach()
        if parent is not None and parent is not boundary_node and len(parent.children) == 1:
            parent.delete(preserve_branch_length=True)
    return boundary_node


def build_skeleton_tree(t, subtree_taxa_map, n_representatives=2):
    """Reduce each analysed subtree to a few representatives in the parent tree.
    """
    representative_map = {}
    skipped = []
    for subtree_name, taxa in subtree_taxa_map.items():
        taxa = [n for n in taxa if n]
        if not taxa:
            continue
        try:
            mrca = t.get_common_ancestor(taxa) if len(taxa) > 1 else t[taxa[0]]
        except Exception as e:
            skipped.append((subtree_name, str(e)))
            continue
        if mrca.is_leaf:
            representative_map[subtree_name] = [mrca.name]
            continue
        mrca_leaf_names = set(l.name for l in mrca.leaves())
        taxa_set = set(taxa) & mrca_leaf_names
        extra = mrca_leaf_names - taxa_set
        if extra:
            print(f"WARNING: subtree '{subtree_name}' isn't perfectly monophyletic in the parent "
                  f"tree ({len(extra)} non-subtree leaf(ves) fall under its MRCA); those are left "
                  f"untouched in the skeleton rather than risk removing something that belongs "
                  f"to a different part of the tree.")
        target_n = min(n_representatives, len(taxa_set))
        treemmer_prune_to_n(mrca, target_n + len(extra), protected_names=extra)
        representative_map[subtree_name] = [l.name for l in mrca.leaves() if l.name in taxa_set]
    return t, representative_map, skipped


def _prune_taxa_by_name(subtree_root, names_to_remove):
    """Detach the given leaf names from subtree_root (in place) and collapse
    any internal nodes left with a single child, walking up the chain."""
    for name in names_to_remove:
        try:
            leaf = subtree_root[name]
        except Exception:
            continue
        if not leaf.is_leaf:
            continue
        parent = leaf.up
        leaf.detach()
        node = parent
        while node is not None and node is not subtree_root and len(node.children) == 1:
            grandparent = node.up
            node.delete(preserve_branch_length=True)
            node = grandparent


def write_tip_ids(subtree, output_path):
    """Write the bare accessions of a subtree's tips, one per line.
    """
    ids = [leaf.name.split('@')[0] for leaf in subtree.leaves()]
    with open(output_path, 'w') as f:
        f.write("\n".join(ids))
    return ids


def extract_fasta_by_ids(master_fasta_path, id_to_output_paths):
    """extract sequences from one master FASTA to many subtree FASTA files.
    """
    open_handles = {}
    found_ids = set()
    try:
        for paths in id_to_output_paths.values():
            for p in paths:
                if p not in open_handles:
                    open_handles[p] = open(p, 'w')
        current_id = None
        current_lines = []
        write_targets = []

        def flush():
            if current_id is not None and write_targets:
                record = ">" + current_id + "\n" + "\n".join(current_lines) + ("\n" if current_lines else "")
                for path in write_targets:
                    open_handles[path].write(record)

        with open(master_fasta_path) as fasta:
            for line in fasta:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    flush()
                    header_tokens = line[1:].split()
                    header_id = header_tokens[0] if header_tokens else line[1:]
                    current_id = header_id
                    current_lines = []
                    write_targets = id_to_output_paths.get(header_id, [])
                    if write_targets:
                        found_ids.add(header_id)
                else:
                    if write_targets:
                        current_lines.append(line)
            flush()
    finally:
        for h in open_handles.values():
            h.close()
    missing = set(id_to_output_paths) - found_ids
    return found_ids, missing


def check_candidates_against_all_families(t, candidate_nodes, ref_prefix='ref|', top_n=3):
    """relate each candidate clade to the curated references.
    """
    all_refs = [l for l in t.leaves() if l.name.startswith(ref_prefix)]
    families = sorted(set(r.name.split('|')[1] for r in all_refs))
    trusted_regions = {fam: get_trusted_family_region(t, fam)[0] for fam in families}

    results = []
    for node in candidate_nodes:
        node_leafnames = set(l.name for l in node.leaves())
        in_trusted = [
            fam for fam, region in trusted_regions.items()
            if region is not None and node_leafnames.issubset(region)
        ]

        scored = []
        for r in all_refs:
            mrca = t.common_ancestor([node, r])
            scored.append((len(list(mrca.leaves())), r.name))
        scored.sort()

        results.append({
            'node': node,
            'in_trusted_region_of': in_trusted,
            'nearest_refs': scored[:top_n],
        })
    return results


def get_leaf_phylum(label, lineage_dict, asg_map_dict):
    """Phylum of a tip, resolved by NCBI identifier and then by curated mapping."""
    id_part = label.split('@')[0]
    ncbi_id = id_part.split('.')[0]
    if ncbi_id in lineage_dict:
        return lineage_dict[ncbi_id]['phylum']
    asg_pattern = re.compile(r'^(.*?)_+\d+$')
    asg_match = asg_pattern.search(id_part)
    if asg_match:
        asg_id = asg_match.group(1)
        asgard_taxon = asg_map_dict.get(asg_id)
        if asgard_taxon:
            return asgard_taxon['phylum']
    return 'Unknown'

def to_nexus(t, output_path):
    """Write a tree to NEXUS with per-tip colours and FigTree display directives.
    """
    nhx = t.write(props=['color'], format_root_node=True)
    tree_str = re.sub(r'\[&&NHX:color=(#[0-9A-Fa-f]+)\]', r'[&color="\1"]', nhx)
    with open(output_path, 'w') as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN TREES;\n")
        f.write(f"\tTREE 1 = [&R] {tree_str}\n")
        f.write("END;\n\n")
        f.write("BEGIN FIGTREE;\n")
        f.write('\tset tipLabels.colorAttribute="color";\n')
        f.write('\tset tipLabels.displayAttribute="Names";\n')
        f.write('\tset nodeLabels.isShown=true;\n')
        f.write('\tset nodeLabels.displayAttribute="label";\n')
        f.write("END;\n")

def get_leaf_domain_phylum(leaf_name, lineage_dict, name_to_taxid, name_to_phylum):
    """Return (domain, phylum) for a tip label.
    """
    if leaf_name.startswith('ref|'):
        parts = leaf_name.split('|')
        family = parts[1] if len(parts) >= 2 else 'Unknown'
        return 'Reference', family

    taxon_name = leaf_name.split('@')[1] if '@' in leaf_name else leaf_name
    taxon_id = get_tax_id(leaf_name, name_to_taxid) if '@' in leaf_name else None
    domain = 'Unknown'
    phylum = 'Unknown'
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
    if not domain:
        domain = 'Unknown'
    if not phylum or "incertae sedis" in phylum.lower() or "unclassified" in phylum.lower():
        phylum = "Unknown"
    return domain, phylum


def _cumulative_top_values(values, thresholds=(75, 90, 100)):
    """Smallest set of values accounting for at least each given percentage, to identify 
    the predominant sister and outgroup donor taxonomy.
    """
    total = len(values)
    if total == 0:
        return {f"{t}%": [] for t in thresholds}
    counts = Counter(values)
    ranked = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    result = {}
    for t in thresholds:
        target = (t / 100.0) * total
        cum = 0
        chosen = []
        for val, cnt in ranked:
            if cum >= target:
                break
            chosen.append(val)
            cum += cnt
        result[f"{t}%"] = chosen
    return result


def summarize_clade_composition(leaf_names, lineage_dict, name_to_taxid, asg_map_dict,
                                 thresholds=(75, 90, 100)):
    """Composition of a tip set, summarised at several cumulative thresholds.
    Returns the taxonomy for each sister and outgroup for each threshold,
    the minimal sets of domains and phyla accounting for that share of the tips.
    """
    name_to_phylum = {data['name']: data['phylum'] for data in asg_map_dict.values()}
    domains, phyla = [], []
    for name in leaf_names:
        d, p = get_leaf_domain_phylum(name, lineage_dict, name_to_taxid, name_to_phylum)
        domains.append(d)
        phyla.append(p)
    domain_summary = _cumulative_top_values(domains, thresholds)
    phylum_summary = _cumulative_top_values(phyla, thresholds)
    return {
        pct: {"domain": domain_summary[pct], "phylum": phylum_summary[pct]}
        for pct in domain_summary
    }


def color_and_export_subtree(tree, lineage_dict, name_to_taxid, asg_map_dict, accesion_db,
                              output_dir, prefix, phylum_to_color=None, color_state=None):
    """Annotate the subtree with taxonomy, assign colours by phylum, and export it.
    """
    if phylum_to_color is None:
        phylum_to_color = {
            'Promethearchaeati': '#ffb347',
            'Unknown': '#708090'
        }
    if color_state is None:
        color_state = {'index': 0}
    color_palette = [
        '#4b99e7', '#fa9200', '#236a80', '#b5b144', '#99566e',
        '#74a8ab', '#74ab85', '#a9011d', '#d00842', '#5b4cf2',
        '#574964', '#232f72', '#9a6324', '#bf360c', '#800000',
        "#1f3b1b", '#808000', '#ffd8b1', '#000075', '#808080'
    ]
    name_to_phylum = {data['name']: data['phylum'] for data in asg_map_dict.values()}
    taxa_path = output_dir / f'{prefix}_taxa.txt'
    nexus_path = output_dir / f'{prefix}.nex'
    with open(taxa_path, 'w') as out_file:
        for leaf in tree.leaves():
            leaf_name = leaf.name
            if '@' not in leaf_name and not leaf_name.startswith('ref|'):
                leaf_name = get_species_name(leaf_name, lineage_dict, accesion_db, asg_map_dict)
                leaf.name = leaf_name
            taxon_name = leaf_name.split('@')[1] if '@' in leaf_name else leaf_name
            domain, phylum = get_leaf_domain_phylum(leaf_name, lineage_dict, name_to_taxid, name_to_phylum)
            leaf.name = leaf.name + '{Phylum: ' + phylum + '}'
            if phylum not in phylum_to_color:
                if color_state['index'] < len(color_palette):
                    phylum_to_color[phylum] = color_palette[color_state['index']]
                    color_state['index'] += 1
                else:
                    import random
                    phylum_to_color[phylum] = "#{:06x}".format(random.randint(0, 0xFFFFFF))
            color = phylum_to_color[phylum]
            leaf.add_prop('color', color)
            out_file.write(f"{taxon_name}\t{phylum}\t{domain}\n")
    to_nexus(tree, nexus_path)
    return phylum_to_color


def select_monophyletic_asgard_subtrees(t, focal_node_dict, focal_list, lineage_dict, name_to_taxid,
                                         asg_map_dict, max_subtree_size=450, max_sister_taxa=300,
                                         composition_thresholds=(75, 90, 100),
                                         sister_absorb_threshold=90.0, focal_phylum='Promethearchaeati'):
    """Build, downsample, characterise and rank a context subtree per focal clade.
    Step 1. Enlarge the clade by absorbing sister clades dominated by the focal phylum,
       then discard enlarged clades that are duplicates of, or nested within,
       another, so that each set of focal sequences is analysed once.
    Step 2. Extract the subtree rooted at the grandparent (focal clade, sisters,
       outgroups) and work on a copy, leaving the input tree intact for the
       remaining candidates.
    Step 3. Downsample any sister clade exceeding `max_sister_taxa`, then the whole
       subtree if it exceeds `max_subtree_size`, protecting the focal tips so that
       size reduction never removes the sequences under study. Both the original and
       retained focal counts are recorded, so any loss is visible.
    Step 4. Count focal tips and distinct non-focal phyla, summarise the composition of
       the sister and outgroup contexts, and classify the transfer event.
    Step 5. Rank eligible candidates by focal-clade support, then number of focal
       sequences, then number of non-focal phyla represented - that is, best
       supported and most informative first.
    Step 6. Enforce exclusivity between selected subtrees: context tips already used by a
       higher-ranked subtree are pruned from the current one (focal tips are never
       pruned), and all composition statistics and the event classification are then
       recomputed for the reduced subtree.
    """
    all_candidates = []
    focal_ids = {tax_id for tax_id, name in focal_list}
    expanded_nodes = {
        num: absorb_dominant_sister_clades(
            node, lineage_dict, name_to_taxid, asg_map_dict,
            target_phylum=focal_phylum, min_dominant_pct=sister_absorb_threshold
        )[0]
        for num, node in focal_node_dict.items()
    }
    expanded_ids = set(id(n) for n in expanded_nodes.values())
    seen_ids = set()
    deduped_nodes = []
    for node in expanded_nodes.values():
        if id(node) in seen_ids:
            continue
        if any(id(anc) in expanded_ids for anc in node.ancestors()):
            continue
        seen_ids.add(id(node))
        deduped_nodes.append(node)
    focal_node_dict = {i: node for i, node in enumerate(deduped_nodes)}

    name_to_phylum = {data['name']: data['phylum'] for data in asg_map_dict.values()}

    for num, focal_node in focal_node_dict.items():
        focal_support = getattr(focal_node, 'support', 1.0)
        if focal_support is None:
            focal_support = 1.0
        focal_size = len(list(focal_node.leaves()))
        subtree_root, sister_nodes, outgroup_nodes = get_subtree_with_context(focal_node)
        subtree_copy = subtree_root.copy()

        def _find_by_leafset(root, target_leafnames):
            for n in root.traverse():
                if set(l.name for l in n.leaves()) == target_leafnames:
                    return n
            return None

        sister_copies = [
            _find_by_leafset(subtree_copy, set(l.name for l in s.leaves()))
            for s in sister_nodes
        ]
        sister_copies = [s for s in sister_copies if s is not None]
        outgroup_copies = [
            _find_by_leafset(subtree_copy, set(l.name for l in o.leaves()))
            for o in outgroup_nodes
        ]
        outgroup_copies = [o for o in outgroup_copies if o is not None]

        sister_size_before = sum(len(list(s.leaves())) for s in sister_copies)
        outgroup_size = sum(len(list(o.leaves())) for o in outgroup_copies)
        for s in sister_copies:
            if len(list(s.leaves())) > max_sister_taxa:
                treemmer_prune_to_n(s, max_sister_taxa)

        sister_size_after = sum(len(list(s.leaves())) for s in sister_copies)
        subtree_taxa = [leaf.name for leaf in subtree_copy.leaves()]
        total_size = len(subtree_taxa)
        sister_taxa = [leaf.name for s in sister_copies for leaf in s.leaves()]
        outgroup_taxa = [leaf.name for o in outgroup_copies for leaf in o.leaves()]
        focal_leaf_names = set(l.name for l in focal_node.leaves())
        truncated_for_size = False
        if total_size > max_subtree_size:
            treemmer_prune_to_n(subtree_copy, max_subtree_size, protected_names=focal_leaf_names)
            truncated_for_size = True
            subtree_taxa = [leaf.name for leaf in subtree_copy.leaves()]
            total_size = len(subtree_taxa)
            sister_taxa = [leaf.name for s in sister_copies for leaf in s.leaves()]
            outgroup_taxa = [leaf.name for o in outgroup_copies for leaf in o.leaves()]
            sister_size_after = len(sister_taxa)
            outgroup_size = len(outgroup_taxa)

        focal_size_retained = len(focal_leaf_names & set(subtree_taxa))
        focal_taxa_retained = [n for n in subtree_taxa if n in focal_leaf_names]
        n_focal_removed_for_size = focal_size - focal_size_retained

        eligible = total_size <= max_subtree_size
        asgard_count = phyla_count = 0
        if eligible:
            asgard_count = sum(1 for leaf in subtree_copy.leaves()
                                if leaf.name.split('@')[0].split('.')[0] in focal_ids)
            non_focal_phyla = set()
            for leaf in subtree_copy.leaves():
                ncbi_id = leaf.name.split('@')[0].split('.')[0]
                if ncbi_id not in focal_ids:
                    phylum = get_leaf_phylum(leaf.name, lineage_dict, asg_map_dict)
                    if phylum != 'Unknown':
                        non_focal_phyla.add(phylum)
            phyla_count = len(non_focal_phyla)

        sister_summary = _group_domain_phylum_summary(sister_taxa, lineage_dict, name_to_taxid, name_to_phylum)
        outgroup_summary = _group_domain_phylum_summary(outgroup_taxa, lineage_dict, name_to_taxid, name_to_phylum)
        event_type = classify_transfer_event(
            sister_summary, len(sister_taxa), outgroup_summary, len(outgroup_taxa),
            focal_phylum=focal_phylum,
        )

        candidate = {
            'clade_id': num,
            'focal_node': focal_node,
            'subtree': subtree_copy,
            'taxa': subtree_taxa,
            'focal_taxa': focal_taxa_retained,
            'focal_size': focal_size_retained,
            'focal_size_original': focal_size,
            'truncated_for_size': truncated_for_size,
            'n_focal_removed_for_size': n_focal_removed_for_size,
            'focal_support': focal_support,
            'sister_taxa': sister_taxa,
            'outgroup_taxa': outgroup_taxa,
            'sister_size_before': sister_size_before,
            'sister_size_after': sister_size_after,
            'outgroup_size': outgroup_size,
            'total_size': total_size,
            'pruned': sister_size_before > sister_size_after,
            'eligible': eligible,
            'selected': False,
            'asgard_count': asgard_count,
            'phyla_count': phyla_count,
            'context_empty': (len(sister_taxa) == 0 and len(outgroup_taxa) == 0),
            'sister_dominant_domain': sister_summary['dominant_domain'],
            'sister_dominant_domain_pct': sister_summary['dominant_domain_pct'],
            'sister_dominant_phylum': sister_summary['dominant_phylum'],
            'sister_dominant_phylum_pct': sister_summary['dominant_phylum_pct'],
            'sister_mixed': sister_summary['mixed'],
            'outgroup_dominant_domain': outgroup_summary['dominant_domain'],
            'outgroup_dominant_domain_pct': outgroup_summary['dominant_domain_pct'],
            'outgroup_dominant_phylum': outgroup_summary['dominant_phylum'],
            'outgroup_dominant_phylum_pct': outgroup_summary['dominant_phylum_pct'],
            'outgroup_mixed': outgroup_summary['mixed'],
            'event_type': event_type,
            'composition': {
                'sister': summarize_clade_composition(
                    sister_taxa, lineage_dict, name_to_taxid, asg_map_dict, composition_thresholds
                ),
                'outgroup': summarize_clade_composition(
                    outgroup_taxa, lineage_dict, name_to_taxid, asg_map_dict, composition_thresholds
                ),
            },
        }
        all_candidates.append(candidate)
    eligible_candidates = [c for c in all_candidates if c['eligible']]
    eligible_candidates.sort(key=lambda c: (c['focal_support'], c['asgard_count'], c['phyla_count']), reverse=True)

    selected = []
    used_taxa = set()
    for c in eligible_candidates:
        focal_names = set(leaf.name for leaf in c['focal_node'].leaves())
        overlap = (set(c['taxa']) & used_taxa) - focal_names
        if overlap:
            _prune_taxa_by_name(c['subtree'], overlap)
            c['taxa'] = [n for n in c['taxa'] if n not in overlap]
            c['sister_taxa'] = [n for n in c['sister_taxa'] if n not in overlap]
            c['outgroup_taxa'] = [n for n in c['outgroup_taxa'] if n not in overlap]
            c['total_size'] = len(c['taxa'])
            c['sister_size_after'] = len(c['sister_taxa'])
            c['outgroup_size'] = len(c['outgroup_taxa'])
            c['pruned'] = True
            c['n_taxa_removed_for_exclusivity'] = len(overlap)
            c['context_empty'] = (len(c['sister_taxa']) == 0 and len(c['outgroup_taxa']) == 0)
            c['asgard_count'] = sum(
                1 for name in c['taxa'] if name.split('@')[0].split('.')[0] in focal_ids
            )
            non_focal_phyla = set()
            for name in c['taxa']:
                ncbi_id = name.split('@')[0].split('.')[0]
                if ncbi_id not in focal_ids:
                    phylum = get_leaf_phylum(name, lineage_dict, asg_map_dict)
                    if phylum != 'Unknown':
                        non_focal_phyla.add(phylum)
            c['phyla_count'] = len(non_focal_phyla)
            sister_summary = _group_domain_phylum_summary(c['sister_taxa'], lineage_dict, name_to_taxid, name_to_phylum)
            outgroup_summary = _group_domain_phylum_summary(c['outgroup_taxa'], lineage_dict, name_to_taxid, name_to_phylum)
            c['sister_dominant_domain'] = sister_summary['dominant_domain']
            c['sister_dominant_domain_pct'] = sister_summary['dominant_domain_pct']
            c['sister_dominant_phylum'] = sister_summary['dominant_phylum']
            c['sister_dominant_phylum_pct'] = sister_summary['dominant_phylum_pct']
            c['sister_mixed'] = sister_summary['mixed']
            c['outgroup_dominant_domain'] = outgroup_summary['dominant_domain']
            c['outgroup_dominant_domain_pct'] = outgroup_summary['dominant_domain_pct']
            c['outgroup_dominant_phylum'] = outgroup_summary['dominant_phylum']
            c['outgroup_dominant_phylum_pct'] = outgroup_summary['dominant_phylum_pct']
            c['outgroup_mixed'] = outgroup_summary['mixed']
            c['event_type'] = classify_transfer_event(
                sister_summary, len(c['sister_taxa']), outgroup_summary, len(c['outgroup_taxa']),
                focal_phylum=focal_phylum,
            )
            c['composition'] = {
                'sister': summarize_clade_composition(
                    c['sister_taxa'], lineage_dict, name_to_taxid, asg_map_dict, composition_thresholds
                ),
                'outgroup': summarize_clade_composition(
                    c['outgroup_taxa'], lineage_dict, name_to_taxid, asg_map_dict, composition_thresholds
                ),
            }
        used_taxa.update(c['taxa'])
        selected.append(c)

    for i, c in enumerate(selected):
        c['selected'] = True
        c['subtree_rank'] = i + 1

    return selected, all_candidates


def json_safe_candidate(c):
    """ make the tip entries safe for it to be serialised to JSON."""
    return {k: v for k, v in c.items() if k not in ('focal_node', 'subtree')}


def write_full_report(output_path, run_info, trusted_clusters_report, all_candidates):
    """Write the full run report: parameters, trusted-region clusters
    and every candidate clade, selected or not."""
    report = {
        "run_info": run_info,
        "trusted_family_clusters": trusted_clusters_report,
        "candidate_clades": [json_safe_candidate(c) for c in all_candidates],
    }
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    return report

def _dominant_value(values):
    if not values:
        return 'NA', 0.0
    counts = Counter(values)
    val, cnt = counts.most_common(1)[0]
    return val, round(100.0 * cnt / len(values), 1)


def _group_domain_phylum_summary(leaf_names, lineage_dict, name_to_taxid, name_to_phylum):
    """Dominant domain and phylum of a tip set, with their percentage shares.
    """
    domains, phyla = [], []
    for name in leaf_names:
        d, p = get_leaf_domain_phylum(name, lineage_dict, name_to_taxid, name_to_phylum)
        domains.append(d)
        phyla.append(p)
    if not domains:
        return {
            'dominant_domain': 'NA', 'dominant_domain_pct': 'NA',
            'dominant_phylum': 'NA', 'dominant_phylum_pct': 'NA',
            'mixed': 'NA',
        }
    dom_domain, dom_domain_pct = _dominant_value(domains)
    dom_phylum, dom_phylum_pct = _dominant_value(phyla)
    return {
        'dominant_domain': dom_domain, 'dominant_domain_pct': dom_domain_pct,
        'dominant_phylum': dom_phylum, 'dominant_phylum_pct': dom_phylum_pct,
        'mixed': dom_domain_pct < 100,
    }


def classify_transfer_event(sister_summary, sister_size, outgroup_summary, outgroup_size,
                             focal_phylum='Promethearchaeati', majority_threshold=50.0):
    """Classify the evolutionary event implied by a clade's phylogenetic context.
    - no context at all                      -> Ambiguous (no context)
    - dominant domain unresolved             -> Ambiguous (unresolved)
    - dominant domain below the majority
      threshold                              -> Ambiguous (mixed context)
    - dominant domain is Archaea             -> Vertical Gene Transfer
    - dominant domain is Reference           -> Ambiguous (reference-only context)
    - otherwise                              -> Horizontal Gene Transfer
    """
    if sister_size > 0:
        summary = sister_summary
    elif outgroup_size > 0:
        summary = outgroup_summary
    else:
        return 'Ambiguous (no context)'

    dom_domain = summary['dominant_domain']
    dom_domain_pct = summary['dominant_domain_pct']
    dom_phylum = summary['dominant_phylum']

    if dom_domain in ('NA', 'Unknown') or dom_domain_pct == 'NA':
        return 'Ambiguous (unresolved)'
    if dom_domain_pct < majority_threshold:
        return 'Ambiguous (mixed context)'
    if dom_phylum == focal_phylum:
        return 'Ambiguous (focal phylum)'
    if dom_domain == 'Archaea':
        return 'Vertical Gene Transfer'
    if dom_domain == 'Reference':
        return 'Ambiguous (reference-only context)'
    return 'Horizontal Gene Transfer'


def _n_unique_species(names):
    """Number of distinct source organisms in a set of tip labels, to check for taxonomically broad sister/outgroup
    """
    species = set()
    for n in names:
        species.add(n.split('@', 1)[1] if '@' in n else n)
    return len(species)


def write_candidate_overview_tsv(all_candidates, output_path, lineage_dict, name_to_taxid,
                                  asg_map_dict, family_type='NA'):
    """Write the candidate review table.
    """
    columns = [
        "clade_id", "eligible", "selected", "subtree_rank", "context_empty",
        "family_type", "event_type",
        "truncated_for_size", "focal_size_original", "n_focal_removed_for_size",
        "focal_size", "focal_spp", "focal_support",
        "sister_size", "sister_spp", "sister_dominant_domain", "sister_dominant_domain_pct",
        "sister_dominant_phylum", "sister_dominant_phylum_pct", "sister_mixed",
        "outgroup_size", "outgroup_spp", "outgroup_dominant_domain", "outgroup_dominant_domain_pct",
        "outgroup_dominant_phylum", "outgroup_dominant_phylum_pct", "outgroup_mixed",
        "total_size",
        "focal_seqs", "sister_seqs", "donor_seqs",
    ]

    def _fmt(x, ndigits=4):
        if x is None:
            return "NA"
        if isinstance(x, float):
            return f"{x:.{ndigits}f}"
        return str(x)

    with open(output_path, 'w') as f:
        f.write("\t".join(columns) + "\n")
        for c in sorted(all_candidates, key=lambda x: x['clade_id']):
            focal_taxa = c.get('focal_taxa')
            if focal_taxa is None and c.get('focal_node') is not None:
                focal_taxa = [leaf.name for leaf in c['focal_node'].leaves()]
            focal_taxa = focal_taxa or []
            sister_taxa = c.get('sister_taxa', [])
            outgroup_taxa = c.get('outgroup_taxa', [])
            donor_seqs = c.get('donor_seqs', outgroup_taxa)

            row = [
                c['clade_id'], c['eligible'], c['selected'], c.get('subtree_rank', ''),
                c.get('context_empty', False),
                family_type, c.get('event_type', 'NA'),
                c.get('truncated_for_size', False), c.get('focal_size_original', 'NA'),
                c.get('n_focal_removed_for_size', 0),
                c.get('focal_size', 'NA'), _n_unique_species(focal_taxa), _fmt(c.get('focal_support'), 3),
                c.get('sister_size_after', 'NA'), _n_unique_species(sister_taxa),
                c.get('sister_dominant_domain', 'NA'), _fmt(c.get('sister_dominant_domain_pct'), 1),
                c.get('sister_dominant_phylum', 'NA'), _fmt(c.get('sister_dominant_phylum_pct'), 1),
                c.get('sister_mixed', 'NA'),
                c.get('outgroup_size', 'NA'), _n_unique_species(outgroup_taxa),
                c.get('outgroup_dominant_domain', 'NA'), _fmt(c.get('outgroup_dominant_domain_pct'), 1),
                c.get('outgroup_dominant_phylum', 'NA'), _fmt(c.get('outgroup_dominant_phylum_pct'), 1),
                c.get('outgroup_mixed', 'NA'),
                c.get('total_size', 'NA'),
                ";".join(focal_taxa), ";".join(sister_taxa), ";".join(donor_seqs),
            ]
            f.write("\t".join(_fmt(x) if not isinstance(x, str) else x for x in row) + "\n")