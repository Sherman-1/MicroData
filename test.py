from Bio import SeqIO 



fasta = SeqIO.parse("representative_sequences.fasta","fasta")


lenghts = [ len(record.seq) for record in fasta ] 



import matplotlib.pyplot as plt 


plt.figure()
plt.hist(lenghts) 
plt.xlabel("Sequence length")
plt.ylabel("Frequency")
plt.title("Sequence length distribution")
plt.savefig("sequence_length_distribution.png")
plt.close()