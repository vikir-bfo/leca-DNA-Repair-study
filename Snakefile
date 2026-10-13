TAXON = glob_wildcards("prot_db/{taxon}.fasta").taxon
K_NUM = glob_wildcards("iter3/hmm_profile/{k_num}.fasta").k_num

wildcard_constraints:
    k_num = r"K\d+",
    taxon = r"[^/]+"

rule all:
    input:
        expand("iter4/hmm_profile/{k_num}.fasta", k_num=K_NUM)

rule hmm_build:
    input:
        fasta = "iter3/hmm_profile/{k_num}.fasta"
    output:
        hmm = "iter3/hmm_profile/{k_num}.hmm"
    shell:
        "hmmbuild {output.hmm} {input.fasta}"

rule hmm_cat:
    input:
        knum = expand("iter3/hmm_profile/{k_num}.hmm", k_num=K_NUM)
    output:
        hmm = "iter3/hmm_profile/KEGG_Master.hmm"
    shell:
        "cat {input.knum} > {output.hmm}"

rule hmm_search:
    input:
        hmm = "iter3/hmm_profile/KEGG_Master.hmm",
        fasta = "prot_db/{taxon}.fasta"
    output:
        tbl = "iter3/hmm_profile/{taxon}.tsv"
    threads: 38
    shell:
        "hmmsearch -E 1e-5 --cpu 38 --tblout {output.tbl} {input.hmm} {input.fasta}"

rule hmm_transpose:
    input:
        expand("iter3/hmm_profile/{taxon}.tsv", taxon=TAXON)
    output:
        expand("iter3/hits/{k_num}.tsv", k_num=K_NUM)
    script:
        "scripts/HMM_transposing.py"

rule best_hits:
    input:
        expand("iter3/hits/{k_num}.tsv", k_num=K_NUM)
    output:
        "iter3/hits/filtered_hits.tsv"
    script:
        "scripts/Best_Unique_Hits.py"

rule group_knum:
    input:
        "iter3/hits/filtered_hits.tsv"
    output:
        "iter3/best_hits/{k_num}.tsv"
    script:
        "scripts/group_knum.py"

rule ext_seq:
    input:
        hits = "iter3/best_hits/{k_num}.tsv",
        fasta = "iter3/all_taxa.fasta"
    output:
        fasta = "iter3/best_hits/{k_num}.fasta",
        index = temp("iter3/best_hits/{k_num}_index.txt")
    shell:
        """
        cut -f1 {input.hits} > {output.index}
        seqtk subseq {input.fasta} {output.index} >> {output.fasta}  
        sed -i 's/*//g' {output.fasta}
        """

rule cd_hit:
    input:
        fasta = "iter3/best_hits/{k_num}.fasta"
    output:
        clstr = "iter3/clustered_hits/{k_num}.fasta.clstr",
        fasta = "iter3/clustered_hits/{k_num}.fasta"
    threads: 38
    shell:
        "cd-hit -i {input.fasta} -o {output.fasta} -c 0.8 -n 5 -M 0 -d 0 -T 38"

rule mafft:
    input:
        fasta = "iter3/clustered_hits/{k_num}.fasta"
    output:
        fasta = "iter3/aligned_hits/{k_num}.fasta"
    threads: 38
    shell:
        "mafft --auto --anysymbol --thread 38 {input.fasta} > {output.fasta}"

rule trimal:
    input:
        fasta = "iter3/aligned_hits/{k_num}.fasta"
    output:
        fasta = "iter3/aligned_hits/trimmed_hits/{k_num}_trimmed.fasta"
    shell:
        "trimal -in {input.fasta} -out {output.fasta} -gt 0.1"

rule fasttree:
    input:
        fasta = "iter3/aligned_hits/trimmed_hits/{k_num}_trimmed.fasta"
    output:
        tree = "iter3/phylotree/{k_num}.nwk"
    threads: 38
    shell:
q        "FastTreeMP -mlnni 4 -out {output.tree} {input.fasta}"

rule treeshrink:
    input:
        tree = "iter3/phylotree/{k_num}.nwk"
    output:
        tree = "iter3/phylotree/treeshrunk/{k_num}_shrunk.nwk",
        logs = "iter3/phylotree/treeshrunk/logs/{k_num}.txt"
    shell:
        """
        run_treeshrink.py -t {input.tree} -c -q 0.05 -m per-gene && \
        mv iter3/phylotree/{wildcards.k_num}_treeshrink/{wildcards.k_num}_0.05.nwk {output.tree} && \
        mv iter3/phylotree/{wildcards.k_num}_treeshrink/{wildcards.k_num}_RS_0.05.txt {output.logs}
        """

rule prune_tree:
    input:
        tree = "iter3/phylotree/treeshrunk/{k_num}_shrunk.nwk",
        guide = "iter3/PhyloGuide.tsv"
    output:
        tree = "iter3/phylotree/pruned/{k_num}_pruned.nwk"
    conda: "envs/ete4.yaml"
    script:
        "scripts/PhyloParse.py"
        
rule seq_dict:
    input:
        fasta = "iter3/all_taxa.fasta"
    output:
        pkl = "iter3/all_taxa.pkl"
    script:
        "scripts/build_seq_dict.py"

rule write_fasta:
    input:
        tree = "iter3/phylotree/pruned/{k_num}_pruned.nwk",
        seq_dict = "iter3/all_taxa.pkl"
    output:
        fasta = "iter3/pruned_seq/{k_num}.fasta"
    conda: "envs/ete4.yaml"
    script:
        "scripts/PrunedSeq.py"
        
rule hmm_seq_align:
    input:
        fasta = "iter3/pruned_seq/{k_num}.fasta"
    output:
        fasta = "iter4/hmm_profile/{k_num}.fasta"
    threads: 38    
    shell:
        "mafft --auto --anysymbol --thread 38 {input.fasta} > {output.fasta}"
