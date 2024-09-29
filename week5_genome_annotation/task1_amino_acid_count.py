# Define the amino acid sequence (without the stop codon '*')
amino_acid_sequence = "KVRMFTSELDIMLSVNGPADQIKYFCRHWT"  # 28 amino acids

# Count the number of amino acids (excluding stop codon)
num_amino_acids = len(amino_acid_sequence)

# Calculate the total number of bases (including stop codon)
total_bases = (num_amino_acids * 3) + 3  # +3 for the stop codon

# Print the results
print(f"Number of amino acids: {num_amino_acids}")
print(f"Total number of bases (including stop codon): {total_bases}")

