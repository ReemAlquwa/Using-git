# Define the amino acid sequence
amino_acid_sequence = "KVRMFTSELDIMLSVNGPADQIKYFCRHWT"

# Calculate number of amino acids
num_amino_acids = len(amino_acid_sequence)

# Calculate number of bases (3 bases per amino acid + 3 bases for stop codon)
num_bases = (num_amino_acids + 1) * 3

print(f"Number of amino acids: {num_amino_acids}")
print(f"Number of bases (including stop codon): {num_bases}")

(myenv) [alquwara@login509-02-l ~]$ nano task1_amino_acid_count.py
(myenv) [alquwara@login509-02-l ~]$ (myenv) [alquwara@login509-02-l ~]$ python task1_amino_acid_count.py
Number of amino acids: 30
Total number of bases (including stop codon): 93
