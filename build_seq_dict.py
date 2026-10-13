from Bio import SeqIO
import pickle

prot_seq = {record.id: record for record in SeqIO.parse(snakemake.input.fasta, "fasta")}
with open(snakemake.output.pkl, "wb") as f:
    pickle.dump(prot_seq, f)