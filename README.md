'Tracing the Evolution of DNA repair machinery from Asgard archaea to Eukaryotes' -- analysis scripts
These are scripts that were used in the study. They cover the homology search post processing; The gene-tree analysis and the LECA repertoire
The requirements are python > 3.9; pandas; ete4; biopython; snakemake along with the NCBI taxonomy file fullnamelineage.dmp, 
accession2taxid.tsv made into a SQL database; cleaned_asgard+thermo_mapping.tsv for the mapping MAG taxa, extracted_taxa.tsv
Scripts:
HMM_transposing.py -- Transposes the HMM output from per-genome to per-profile.
Best_Unique_Hits.py -- collapses the genome derived hits to a representative per-gene
Unique_HMMs.py -- Reduces the per profile tables to a best hit per protein sequence table
subtree_selection.py -- library that holds all the functions for the gene tree analysis of asgards
subtree_main.upd.py -- the frontend for getting the selected subtreed for a parent tree
python subtree_main.upd.py \
    --tree        
    --lineage     fullnamelineage.dmp \
    --asg_lineage cleaned_asgard+thermo_mapping.tsv \
    --rep_db      rep_dict.db \
    --out         
colour_subtree.order.py -- Rewrites tip labels of one selected subtree with organism and taxonomy, exports a NEXUS file that is coloured.
python colour_subtree.order.py \
    --tree   
    --output 
domain_mapping.py -- draws subtrees with the InterPro domains coloured by accession
python domain_mapping.py \
    --tree    
    --tsv     interproscan_output.tsv \
    --db      Pfam \
    --outfile 
tree_prune.py -- Removes the out-paralogues and in-paralogues from a eukaryotic gene tree
