from Bio import SeqIO
import pickle
""" To build a pickle dump with all the fasta files in them"""
prot_seq = {record.id: record for record in SeqIO.parse(snakemake.input.fasta, "fasta")}
with open(snakemake.output.pkl, "wb") as f:
    pickle.dump(prot_seq, f)
