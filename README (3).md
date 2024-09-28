
# Python tool to find genes in a genome using the FASTA format

 A Python script that takes a FASTA file as input and outputs regions between a start codon (ATG) and stop codons (TAA, TAG, TGA) in three reading frames.## Set Up a Python Virtual Environment  

# Navigate to your gene_finder directory:

cd ~/gene_finder

# Create a virtual environment using Python's venv module:
python3 -m venv myenv


# Activate the virtual environment:
source myenv/bin/activate


# Install Biopython Inside the Virtual Environment
pip install biopython



## Log in to Ibex

#  gene_finder directory create it:

ssh Alquwara@ilogin.ibex.kaust.edu.sa



mkdir ~/gene_finder
cd ~/gene_finder

cd ~/gene_finder/Reem/



#  Activating the Python Environment
source ~/myenv/bin/activate

## Task 1: Setup and Run - Open the nano Editor
 # To create or edit the Python script

cd ~/gene_finder

nano gene_finder_task1.py



## Writing the  Task 1 code Python in (nano):

import sys

def read_fasta(filename):
    """Reads a single FASTA file and returns the sequence."""
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

def find_genes(sequence):
    """Finds genes in the forward reading frames."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    genes = []

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        gene = sequence[i:j+3]
                        genes.append(gene)
                        break
    return genes

def filter_genes_by_length(genes, min_length=100):
    """Filter genes by length (minimum codon length)."""
    return [gene for gene in genes if len(gene) >= min_length * 3]

def main():
    if len(sys.argv) != 2:
        print("Usage: python gene_finder_task1.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

    # Find genes in forward reading frames
    genes = find_genes(sequence)
    filtered_genes = filter_genes_by_length(genes)
    output_file = fasta_file.replace('.fna', '_task1_output.txt')

    # Save the filtered genes to the output file
    with open(output_file, 'w') as f:
        for gene in filtered_genes:
            f.write(f"Filtered Gene: {gene}\n")

if __name__ == "__main__":
    main()



## Save and Exit nano

## To run Task 1:
 python gene_finder_task1.py Reem/GCA_000006745.1_ASM674v1_genomic.fna



# to Show it : 
less Reem/GCA_000006745.1_ASM674v1_genomic_task1_output.txt


# output:

Filtered Gene: ATGATGTGCAGCGGCATCCCATCAATATGGATATGCTCACGCAGAACATCACGGGTGGTACCGGCAATGTCGGTAACGATGGCAGACTCTTTACCTGAAAGCGCATTGAGTAGGCTCGATTTACCCGCATTAGGACGCCCAGCAATCACCACCTTCATCCCTTCGCGCATAATGGCGCCTTGGTTGGCTTCACGGCGCACTGCGGCAAGATTATCTATGATGGTTTGCAGATCAGCGGAAACCTTACCATCGGCCAGAAAATCGATCTCTTCTTCTGGGAAATCAATTGCGGCTTCAACATAG
Filtered Gene: ATGTGCAGCGGCATCCCATCAATATGGATATGCTCACGCAGAACATCACGGGTGGTACCGGCAATGTCGGTAACGATGGCAGACTCTTTACCTGAAAGCGCATTGAGTAGGCTCGATTTACCCGCATTAGGACGCCCAGCAATCACCACCTTCATCCCTTCGCGCATAATGGCGCCTTGGTTGGCTTCACGGCGCACTGCGGCAAGATTATCTATGATGGTTTGCAGATCAGCGGAAACCTTACCATCGGCCAGAAAATCGATCTCTTCTTCTGGGAAATCAATTGCGGCTTCAACATAG
Filtered Gene: ATGCCGCTTTCCCGGTTAATGATCCAACAGTTTCGTAATATTAAAGCCTGTGATATTCGCCTATCAGCAGGCTTTAACTTTCTTATCGGGCCGAATGGCAGCGGTAAAACCAGCGTGCTTGAAGCGATTTATCTGCTCGGCCATGGCCGTTCATTTAAAAGCTCGCTGACCGGACGGATCATTCAAAATGAGTGCTCAGAACTGTTTGTGCATGGCCGGATTTGTGAGCATTCTTTGAGCTCTGATCAATTTGAGCTACCCGTTGGCATTAATAAGCAGCGTGACGGCTCAACTGAGGTTAAAATAGGCGGCCAAACAGGGCAAAAATTGGCGCAATTAGCGCAGATTTTGCCACTGCAATTGATTCATCCGGAAGGCTTTGAACTGCTGACCGATGGACCGAAACAGCGTCGGGCTTTTATCGATTGGGGCGTGTTTCATACTGAGCCCGCTTTTTTCGATGCATGGGGGCGGTTTAAGCGGCTCAGCAAGCAGCGTAATGCACTGCTAAAAAGCGCGCAAAGCTATCGTGAACTCAGTTATTGGGATCAAGAATTAGCCCGTTTGGCTGAGCAGATTGACCAGTGGCGTGAAAGTTACGTCAATCAATTGAAAAATGTGGCAGAGCAGTTATGTCGCACATTTTTGCCAGAATTCGATATCGACTTAAAGTATTATCGAGGTTGGGAAAAAGATCAGCCTTATCAATCGATTCTGGAAAAAAACTTCGAACGGGATCAGCAGTTGGGCTATACCTTTAGCGGGCCTAACAAAGCGGATTTGCGGATTAAAGTGAACGCAACCCCGGTGGAAGATGTGTTGTCGCGGGGACAACTGAAGTTGATGGTGTGTGCGTTGCGCGTGGCGCAAGGGCAGCACCTGACCGAGTTGACGGGAAAGCAATGCATTTATCTTATCGATGATTTTGCTTCCGAATTGGATAGCCTACGTCGCCAACGTTTGGCAGATAGCCTAAAAGGGACGGGGGCGCAGGTTTTTGTAAGTTCTATTACCGAAAGCCAAGTGGCGGACATGTTAGATGAGTCCAGTAAAACATTTCACGTTGCCCACGGTGTAATAGAGCAAGGATAA
Filtered Gene: ATGATCCAACAGTTTCGTAATATTAAAGCCTGTGATATTCGCCTATCAGCAGGCTTTAACTTTCTTATCGGGCCGAATGGCAGCGGTAAAACCAGCGTGCTTGAAGCGATTTATCTGCTCGGCCATGGCCGTTCATTTAAAAGCTCGCTGACCGGACGGATCATTCAAAATGAGTGCTCAGAACTGTTTGTGCATGGCCGGATTTGTGAGCATTCTTTGAGCTCTGATCAATTTGAGCTACCCGTTGGCATTAATAAGCAGCGTGACGGCTCAACTGAGGTTAAAATAGGCGGCCAAACAGGGCAAAAATTGGCGCAATTAGCGCAGATTTTGCCACTGCAATTGATTCATCCGGAAGGCTTTGAACTGCTGACCGATGGACCGAAACAGCGTCGGGCTTTTATCGATTGGGGCGTGTTTCATACTGAGCCCGCTTTTTTCGATGCATGGGGGCGGTTTAAGCGGCTCAGCAAGCAGCGTAATGCACTGCTAAAAAGCGCGCAAAGCTATCGTGAACTCAGTTATTGGGATCAAGAATTAGCCCGTTTGGCTGAGCAGATTGACCAGTGGCGTGAAAGTTACGTCAATCAATTGAAAAATGTGGCAGAGCAGTTATGTCGCACATTTTTGCCAGAATTCGATATCGACTTAAAGTATTATCGAGGTTGGGAAAAAGATCAGCCTTATCAATCGATTCTGGAAAAAAACTTCGAACGGGATCAGCAGTTGGGCTATACCTTTAGCGGGCCTAACAAAGCGGATTTGCGGATTAAAGTGAACGCAACCCCGGTGGAAGATGTGTTGTCGCGGGGACAACTGAAGTTGATGGTGTGTGCGTTGCGCGTGGCGCAAGGGCAGCACCTGACCGAGTTGACGGGAAAGCAATGCATTTATCTTATCGATGATTTTGCTTCCGAATTGGATAGCCTACGTCGCCAACGTTTGGCAGATAGCCTAAAAGGGACGGGGGCGCAGGTTTTTGTAAGTTCTATTACCGAAAGCCAAGTGGCGGACATGTTAGATGAGTCCAGTAAAACATTTCACGTTGCCCACGGTGTAATAGAGCAAGGATAA
## Task 2 -Reverse 
cd ~/gene_finder

nano gene_finder_task2.py


## Writing the  Task 2 code Python in (nano):


import sys

def flipComp(Text):
#Compute the reverse complement of a DNA sequence.
    result = ""
    for letter in Text:
        if letter == 'A':
            result += 'T'
        elif letter == 'T':
            result += 'A'
        elif letter == 'G':
            result += 'C'
        else:
            result += 'G'
    flipped = result[::-1]  # Reverse the sequence
    return flipped

def read_fasta(filename):
#Reads a single FASTA file and returns the sequence.
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

def find_genes(sequence, include_reverse=False):
    """Finds genes in reading frames (with or without reverse complements)."""
    forward_start_codon = 'ATG'
    forward_stop_codons = ['TAA', 'TAG', 'TGA']

    reverse_start_codon = 'CAT'  # Reverse complement of ATG
    reverse_stop_codons = ['TTA', 'CTA', 'TCA']  # Reverse complements of stop codons

    genes = []

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == forward_start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in forward_stop_codons:
                        gene = sequence[i:j+3]
                        genes.append(f"Forward: {gene}")
                        break

    # Check reverse reading frames if required
    if include_reverse:
        rev_sequence = flipComp(sequence)  # Get reverse complement using flipComp
        print(f"Reverse Complement Sequence: {rev_sequence[:50]}...")  # Print first 50 bases of reverse complement

        for frame in range(3):
            for i in range(frame, len(rev_sequence), 3):
                codon = rev_sequence[i:i+3]
                if codon == reverse_start_codon:  # Check for CAT in reverse complement
                    print(f"Found Reverse Start Codon (CAT) at position {i} in frame {frame}")  # Debug statement
                    for j in range(i+3, len(rev_sequence), 3):
                        stop_codon = rev_sequence[j:j+3]
                        if stop_codon in reverse_stop_codons:  # Look for TTA, CTA, TCA
                            gene = rev_sequence[i:j+3]
                            genes.append(f"Reverse: {gene}")
                            print(f"Found Reverse Gene from {i} to {j}")  # Debug statement
                            break

    return genes

def filter_genes_by_length(genes, min_length=100):
    """Filter genes by length (minimum codon length)."""
    filtered_genes = []
    for gene in genes:
        gene_sequence = gene.split(': ')[1]  # Extract the gene sequence
        if len(gene_sequence) >= min_length * 3:  # Apply length filtering
            filtered_genes.append(gene)
    return filtered_genes

def main():
    if len(sys.argv) != 2:
        print("Usage: python gene_finder_task2.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

# Find genes in forward and reverse reading frames
    genes = find_genes(sequence, include_reverse=True)

# Apply filtering by gene length
    filtered_genes = filter_genes_by_length(genes)

    output_file = fasta_file.replace('.fna', '_task3_forward_reverse_output.txt')

# Save both forward and reverse genes to the output file
    with open(output_file, 'w') as f:
        for gene in filtered_genes:
            if "Forward" in gene:
                f.write(f"Filtered Gene (Forward): {gene.split(': ')[1]}\n")
            elif "Reverse" in gene:
                f.write(f"Filtered Gene (Reverse): {gene.split(': ')[1]}\n")

    print(f"Total genes (forward and reverse) written: {len(filtered_genes)}")  # Debug statement

if __name__ == "__main__":
    main()



## Run Task 2 with the following command:

python gene_finder_task2.py Reem/GCA_000006745.1_ASM674v1_genomic.fna

# To show results 

less Reem/GCA_000006745.1_ASM674v1_genomic_task3_forward_reverse_output.txt 


## output
Filtered Gene (Reverse): CATGCGAGCAAAAGTGACACGCATGTCTTGCGCCGTTTTGAGTGGATCTGGGTTTCCATCCACACCTTCTGGGTTGACGTAGATAAGCCCCATCATGACGGCCGCAAGTGGGTTTTCGAGATCGCGCTGACCAGAGTAACGGCTATTTTCACCGCCAGATTTGGCAAGCCACTCTTTTTCAGAACCCCAGTAAATGTCTTTTTCTGGGTGCCAGATGTCTTCGCGGCCAAAGGCAAAACCAAAGGTTTTGAGGCCCATCGATTCATAAGCCATGTTCCCCGCGAGGATCATCAAATCCGCCCAGCTGATTTTGTTACCGTACTTTTGTTTGATAGGCCACAGCAAACGGCGTGCTTTATCTAAGTTGGCGTTATCAGGCCATGAGTTCAGTGGGGCAAAACGTTGGTTACCGGTACCGCCGCCGCCACGACCATCGGCAATACGGTAGGTACCCGCAGAGTGCCATGCCATACGAATCATCAAACCGCCGTAATGTCCCCAGTCGGCCGGCCACCACTCTTGGCTATTGGTCATTAACGCTTTCAGATCGCGTTTTAGCGCTTCTACATCAAGCTTTTTCAGCTCTTCGCGATAGTTGAAATCAGCACCAAGTGGATTGGTTTTGCTGTCATGCTGATGAAGGATGTCTAAGTTTAGCGCTTTAGGCCACCAGTCCATGTTGGACATGCTCGCCGAAGTCAAGCCACCGTGCATGACAGGGCATTGACCGCTTGAACCCGCTTTATTGTGCTCCATTGATTA
Filtered Gene (Reverse): CATGTCTTGCGCCGTTTTGAGTGGATCTGGGTTTCCATCCACACCTTCTGGGTTGACGTAGATAAGCCCCATCATGACGGCCGCAAGTGGGTTTTCGAGATCGCGCTGACCAGAGTAACGGCTATTTTCACCGCCAGATTTGGCAAGCCACTCTTTTTCAGAACCCCAGTAAATGTCTTTTTCTGGGTGCCAGATGTCTTCGCGGCCAAAGGCAAAACCAAAGGTTTTGAGGCCCATCGATTCATAAGCCATGTTCCCCGCGAGGATCATCAAATCCGCCCAGCTGATTTTGTTACCGTACTTTTGTTTGATAGGCCACAGCAAACGGCGTGCTTTATCTAAGTTGGCGTTATCAGGCCATGAGTTCAGTGGGGCAAAACGTTGGTTACCGGTACCGCCGCCGCCACGACCATCGGCAATACGGTAGGTACCCGCAGAGTGCCATGCCATACGAATCATCAAACCGCCGTAATGTCCCCAGTCGGCCGGCCACCACTCTTGGCTATTGGTCATTAACGCTTTCAGATCGCGTTTTAGCGCTTCTACATCAAGCTTTTTCAGCTCTTCGCGATAGTTGAAATCAGCACCAAGTGGATTGGTTTTGCTGTCATGCTGATGAAGGATGTCTAAGTTTAGCGCTTTAGGCCACCAGTCCATGTTGGACATGCTCGCCGAAGTCAAGCCACCGTGCATGACAGGGCATTGACCGCTTGAACCCGCTTTATTGTGCTCCATTGATTA
Filtered Gene (Reverse): CATCATGACGGCCGCAAGTGGGTTTTCGAGATCGCGCTGACCAGAGTAACGGCTATTTTCACCGCCAGATTTGGCAAGCCACTCTTTTTCAGAACCCCAGTAAATGTCTTTTTCTGGGTGCCAGATGTCTTCGCGGCCAAAGGCAAAACCAAAGGTTTTGAGGCCCATCGATTCATAAGCCATGTTCCCCGCGAGGATCATCAAATCCGCCCAGCTGATTTTGTTACCGTACTTTTGTTTGATAGGCCACAGCAAACGGCGTGCTTTATCTAAGTTGGCGTTATCAGGCCATGAGTTCAGTGGGGCAAAACGTTGGTTACCGGTACCGCCGCCGCCACGACCATCGGCAATACGGTAGGTACCCGCAGAGTGCCATGCCATACGAATCATCAAACCGCCGTAATGTCCCCAGTCGGCCGGCCACCACTCTTGGCTATTGGTCATTAACGCTTTCAGATCGCGTTTTAGCGCTTCTACATCAAGCTTTTTCAGCTCTTCGCGATAGTTGAAATCAGCACCAAGTGGATTGGTTTTGCTGTCATGCTGATGAAGGATGTCTAAGTTTAGCGCTTTAGGCCACCAGTCCATGTTGGACATGCTCGCCGAAGTCAAGCCACCGTGCATGACAGGGCATTGACCGCTTGAACCCGCTTTATTGTGCTCCATTGATTA
... 
Filtered Gene (Forward): ATGTGCAGCGGCATCCCATCAATATGGATATGCTCACGCAGAACATCACGGGTGGTACCGGCAATGTCGGTAACGATGGCAGACTCTTTACCTGAAAGCGCATTGAGTAGGCTCGATTTACCCGCATTAGGACGCCCAGCAATCACCACCTTCATCCCTTCGCGCATAATGGCGCCTTGGTTGGCTTCACGGCGCACTGCGGCAAGATTATCTATGATGGTTTGCAGATCAGCGGAAACCTTACCATCGGCCAGAAAATCGATCTCTTCTTCTGGGAAATCAATTGCGGCTTCAACATAG
Filtered Gene (Forward): ATGCCGCTTTCCCGGTTAATGATCCAACAGTTTCGTAATATTAAAGCCTGTGATATTCGCCTATCAGCAGGCTTTAACTTTCTTATCGGGCCGAATGGCAGCGGTAAAACCAGCGTGCTTGAAGCGATTTATCTGCTCGGCCATGGCCGTTCATTTAAAAGCTCGCTGACCGGACGGATCATTCAAAATGAGTGCTCAGAACTGTTTGTGCATGGCCGGATTTGTGAGCATTCTTTGAGCTCTGATCAATTTGAGCTACCCGTTGGCATTAATAAGCAGCGTGACGGCTCAACTGAGGTTAAAATAGGCGGCCAAACAGGGCAAAAATTGGCGCAATTAGCGCAGATTTTGCCACTGCAATTGATTCATCCGGAAGGCTTTGAACTGCTGACCGATGGACCGAAACAGCGTCGGGCTTTTATCGATTGGGGCGTGTTTCATACTGAGCCCGCTTTTTTCGATGCATGGGGGCGGTTTAAGCGGCTCAGCAAGCAGCGTAATGCACTGCTAAAAAGCGCGCAAAGCTATCGTGAACTCAGTTATTGGGATCAAGAATTAGCCCGTTTGGCTGAGCAGATTGACCAGTGGCGTGAAAGTTACGTCAATCAATTGAAAAATGTGGCAGAGCAGTTATGTCGCACATTTTTGCCAGAATTCGATATCGACTTAAAGTATTATCGAGGTTGGGAAAAAGATCAGCCTTATCAATCGATTCTGGAAAAAAACTTCGAACGGGATCAGCAGTTGGGCTATACCTTTAGCGGGCCTAACAAAGCGGATTTGCGGATTAAAGTGAACGCAACCCCGGTGGAAGATGTGTTGTCGCGGGGACAACTGAAGTTGATGGTGTGTGCGTTGCGCGTGGCGCAAGGGCAGCACCTGACCGAGTTGACGGGAAAGCAATGCATTTATCTTATCGATGATTTTGCTTCCGAATTGGATAGCCTACGTCGCCAACGTTTGGCAGATAGCCTAAAAGGGACGGGGGCGCAGGTTTTTGTAAGTTCTATTACCGAAAGCCAAGTGGCGGACATGTTAGATGAGTCCAGTAAAACATTTCACGTTGCCCACGGTGTAATAGAGCAAGGATAA
Filtered Gene (Forward): ATGATCCAACAGTTTCGTAATATTAAAGCCTGTGATATTCGCCTATCAGCAGGCTTTAACTTTCTTATCGGGCCGAATGGCAGCGGTAAAACCAGCGTGCTTGAAGCGATTTATCTGCTCGGCCATGGCCGTTCATTTAAAAGCTCGCTGACCGGACGGATCATTCAAAATGAGTGCTCAGAACTGTTTGTGCATGGCCGGATTTGTGAGCATTCTTTGAGCTCTGATCAATTTGAGCTACCCGTTGGCATTAATAAGCAGCGTGACGGCTCAACTGAGGTTAAAATAGGCGGCCAAACAGGGCAAAAATTGGCGCAATTAGCGCAGATTTTGCCACTGCAATTGATTCATCCGGAAGGCTTTGAACTGCTGACCGATGGACCGAAACAGCGTCGGGCTTTTATCGATTGGGGCGTGTTTCATACTGAGCCCGCTTTTTTCGATGCATGGGGGCGGTTTAAGCGGCTCAGCAAGCAGCGTAATGCACTGCTAAAAAGCGCGCAAAGCTATCGTGAACTCAGTTATTGGGATCAAGAATTAGCCCGTTTGGCTGAGCAGATTGACCAGTGGCGTGAAAGTTACGTCAATCAATTGAAAAATGTGGCAGAGCAGTTATGTCGCACATTTTTGCCAGAATTCGATATCGACTTAAAGTATTATCGAGGTTGGGAAAAAGATCAGCCTTATCAATCGATTCTGGAAAAAAACTTCGAACGGGATCAGCAGTTGGGCTATACCTTTAGCGGGCCTAACAAAGCGGATTTGCGGATTAAAGTGAACGCAACCCCGGTGGAAGATGTGTTGTCGCGGGGACAACTGAAGTTGATGGTGTGTGCGTTGCGCGTGGCGCAAGGGCAGCACCTGACCGAGTTGACGGGAAAGCAATGCATTTATCTTATCGATGATTTTGCTTCCGAATTGGATAGCCTACGTCGCCAACGTTTGGCAGATAGCCTAAAAGGGACGGGGGCGCAGGTTTTTGTAAGTTCTATTACCGAAAGCCAAGTGGCGGACATGTTAGATGAGTCCAGTAAAACATTTCACGTTGCCCACGGTGTAATAGAGCAAGGATAA






##  Pushing the Output Files to GitHub

cd ~/gene_finder

git status




git add Reem/GCA_000006745.1_ASM674v1_genomic_task1_output.txt

git add Reem/GCA_000006745.1_ASM674v1_genomic_task3_forward_reverse_output.txt


# Commit changes 

git commit -m "Added Task 1 and Task 2 output files"
 
 




##  Pushing the Output Files to GitHub

cd ~/gene_finder

git status




git add Reem/GCA_000006745.1_ASM674v1_genomic_task1_output.txt

git add Reem/GCA_000006745.1_ASM674v1_genomic_task3_forward_reverse_output.txt


# Commit changes 

git commit -m "Added Task 1 and Task 2 output files"
 
 
## Task 2 -Reverse 
cd ~/gene_finder

nano gene_finder_task2.py


## Writing the  Task 2 code Python in (nano):


import sys

def flipComp(Text):
#Compute the reverse complement of a DNA sequence.
    result = ""
    for letter in Text:
        if letter == 'A':
            result += 'T'
        elif letter == 'T':
            result += 'A'
        elif letter == 'G':
            result += 'C'
        else:
            result += 'G'
    flipped = result[::-1]  # Reverse the sequence
    return flipped

def read_fasta(filename):
#Reads a single FASTA file and returns the sequence.
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

def find_genes(sequence, include_reverse=False):
    """Finds genes in reading frames (with or without reverse complements)."""
    forward_start_codon = 'ATG'
    forward_stop_codons = ['TAA', 'TAG', 'TGA']

    reverse_start_codon = 'CAT'  # Reverse complement of ATG
    reverse_stop_codons = ['TTA', 'CTA', 'TCA']  # Reverse complements of stop codons

    genes = []

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == forward_start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in forward_stop_codons:
                        gene = sequence[i:j+3]
                        genes.append(f"Forward: {gene}")
                        break

    # Check reverse reading frames if required
    if include_reverse:
        rev_sequence = flipComp(sequence)  # Get reverse complement using flipComp
        print(f"Reverse Complement Sequence: {rev_sequence[:50]}...")  # Print first 50 bases of reverse complement

        for frame in range(3):
            for i in range(frame, len(rev_sequence), 3):
                codon = rev_sequence[i:i+3]
                if codon == reverse_start_codon:  # Check for CAT in reverse complement
                    print(f"Found Reverse Start Codon (CAT) at position {i} in frame {frame}")  # Debug statement
                    for j in range(i+3, len(rev_sequence), 3):
                        stop_codon = rev_sequence[j:j+3]
                        if stop_codon in reverse_stop_codons:  # Look for TTA, CTA, TCA
                            gene = rev_sequence[i:j+3]
                            genes.append(f"Reverse: {gene}")
                            print(f"Found Reverse Gene from {i} to {j}")  # Debug statement
                            break

    return genes

def filter_genes_by_length(genes, min_length=100):
    """Filter genes by length (minimum codon length)."""
    filtered_genes = []
    for gene in genes:
        gene_sequence = gene.split(': ')[1]  # Extract the gene sequence
        if len(gene_sequence) >= min_length * 3:  # Apply length filtering
            filtered_genes.append(gene)
    return filtered_genes

def main():
    if len(sys.argv) != 2:
        print("Usage: python gene_finder_task2.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

# Find genes in forward and reverse reading frames
    genes = find_genes(sequence, include_reverse=True)

# Apply filtering by gene length
    filtered_genes = filter_genes_by_length(genes)

    output_file = fasta_file.replace('.fna', '_task3_forward_reverse_output.txt')

# Save both forward and reverse genes to the output file
    with open(output_file, 'w') as f:
        for gene in filtered_genes:
            if "Forward" in gene:
                f.write(f"Filtered Gene (Forward): {gene.split(': ')[1]}\n")
            elif "Reverse" in gene:
                f.write(f"Filtered Gene (Reverse): {gene.split(': ')[1]}\n")

    print(f"Total genes (forward and reverse) written: {len(filtered_genes)}")  # Debug statement

if __name__ == "__main__":
    main()



## Run Task 2 with the following command:

python gene_finder_task2.py Reem/GCA_000006745.1_ASM674v1_genomic.fna

# To show results 

less Reem/GCA_000006745.1_ASM674v1_genomic_task3_forward_reverse_output.txt 


## output
Filtered Gene (Reverse): CATGCGAGCAAAAGTGACACGCATGTCTTGCGCCGTTTTGAGTGGATCTGGGTTTCCATCCACACCTTCTGGGTTGACGTAGATAAGCCCCATCATGACGGCCGCAAGTGGGTTTTCGAGATCGCGCTGACCAGAGTAACGGCTATTTTCACCGCCAGATTTGGCAAGCCACTCTTTTTCAGAACCCCAGTAAATGTCTTTTTCTGGGTGCCAGATGTCTTCGCGGCCAAAGGCAAAACCAAAGGTTTTGAGGCCCATCGATTCATAAGCCATGTTCCCCGCGAGGATCATCAAATCCGCCCAGCTGATTTTGTTACCGTACTTTTGTTTGATAGGCCACAGCAAACGGCGTGCTTTATCTAAGTTGGCGTTATCAGGCCATGAGTTCAGTGGGGCAAAACGTTGGTTACCGGTACCGCCGCCGCCACGACCATCGGCAATACGGTAGGTACCCGCAGAGTGCCATGCCATACGAATCATCAAACCGCCGTAATGTCCCCAGTCGGCCGGCCACCACTCTTGGCTATTGGTCATTAACGCTTTCAGATCGCGTTTTAGCGCTTCTACATCAAGCTTTTTCAGCTCTTCGCGATAGTTGAAATCAGCACCAAGTGGATTGGTTTTGCTGTCATGCTGATGAAGGATGTCTAAGTTTAGCGCTTTAGGCCACCAGTCCATGTTGGACATGCTCGCCGAAGTCAAGCCACCGTGCATGACAGGGCATTGACCGCTTGAACCCGCTTTATTGTGCTCCATTGATTA
Filtered Gene (Reverse): CATGTCTTGCGCCGTTTTGAGTGGATCTGGGTTTCCATCCACACCTTCTGGGTTGACGTAGATAAGCCCCATCATGACGGCCGCAAGTGGGTTTTCGAGATCGCGCTGACCAGAGTAACGGCTATTTTCACCGCCAGATTTGGCAAGCCACTCTTTTTCAGAACCCCAGTAAATGTCTTTTTCTGGGTGCCAGATGTCTTCGCGGCCAAAGGCAAAACCAAAGGTTTTGAGGCCCATCGATTCATAAGCCATGTTCCCCGCGAGGATCATCAAATCCGCCCAGCTGATTTTGTTACCGTACTTTTGTTTGATAGGCCACAGCAAACGGCGTGCTTTATCTAAGTTGGCGTTATCAGGCCATGAGTTCAGTGGGGCAAAACGTTGGTTACCGGTACCGCCGCCGCCACGACCATCGGCAATACGGTAGGTACCCGCAGAGTGCCATGCCATACGAATCATCAAACCGCCGTAATGTCCCCAGTCGGCCGGCCACCACTCTTGGCTATTGGTCATTAACGCTTTCAGATCGCGTTTTAGCGCTTCTACATCAAGCTTTTTCAGCTCTTCGCGATAGTTGAAATCAGCACCAAGTGGATTGGTTTTGCTGTCATGCTGATGAAGGATGTCTAAGTTTAGCGCTTTAGGCCACCAGTCCATGTTGGACATGCTCGCCGAAGTCAAGCCACCGTGCATGACAGGGCATTGACCGCTTGAACCCGCTTTATTGTGCTCCATTGATTA
Filtered Gene (Reverse): CATCATGACGGCCGCAAGTGGGTTTTCGAGATCGCGCTGACCAGAGTAACGGCTATTTTCACCGCCAGATTTGGCAAGCCACTCTTTTTCAGAACCCCAGTAAATGTCTTTTTCTGGGTGCCAGATGTCTTCGCGGCCAAAGGCAAAACCAAAGGTTTTGAGGCCCATCGATTCATAAGCCATGTTCCCCGCGAGGATCATCAAATCCGCCCAGCTGATTTTGTTACCGTACTTTTGTTTGATAGGCCACAGCAAACGGCGTGCTTTATCTAAGTTGGCGTTATCAGGCCATGAGTTCAGTGGGGCAAAACGTTGGTTACCGGTACCGCCGCCGCCACGACCATCGGCAATACGGTAGGTACCCGCAGAGTGCCATGCCATACGAATCATCAAACCGCCGTAATGTCCCCAGTCGGCCGGCCACCACTCTTGGCTATTGGTCATTAACGCTTTCAGATCGCGTTTTAGCGCTTCTACATCAAGCTTTTTCAGCTCTTCGCGATAGTTGAAATCAGCACCAAGTGGATTGGTTTTGCTGTCATGCTGATGAAGGATGTCTAAGTTTAGCGCTTTAGGCCACCAGTCCATGTTGGACATGCTCGCCGAAGTCAAGCCACCGTGCATGACAGGGCATTGACCGCTTGAACCCGCTTTATTGTGCTCCATTGATTA
... 
Filtered Gene (Forward): ATGTGCAGCGGCATCCCATCAATATGGATATGCTCACGCAGAACATCACGGGTGGTACCGGCAATGTCGGTAACGATGGCAGACTCTTTACCTGAAAGCGCATTGAGTAGGCTCGATTTACCCGCATTAGGACGCCCAGCAATCACCACCTTCATCCCTTCGCGCATAATGGCGCCTTGGTTGGCTTCACGGCGCACTGCGGCAAGATTATCTATGATGGTTTGCAGATCAGCGGAAACCTTACCATCGGCCAGAAAATCGATCTCTTCTTCTGGGAAATCAATTGCGGCTTCAACATAG
Filtered Gene (Forward): ATGCCGCTTTCCCGGTTAATGATCCAACAGTTTCGTAATATTAAAGCCTGTGATATTCGCCTATCAGCAGGCTTTAACTTTCTTATCGGGCCGAATGGCAGCGGTAAAACCAGCGTGCTTGAAGCGATTTATCTGCTCGGCCATGGCCGTTCATTTAAAAGCTCGCTGACCGGACGGATCATTCAAAATGAGTGCTCAGAACTGTTTGTGCATGGCCGGATTTGTGAGCATTCTTTGAGCTCTGATCAATTTGAGCTACCCGTTGGCATTAATAAGCAGCGTGACGGCTCAACTGAGGTTAAAATAGGCGGCCAAACAGGGCAAAAATTGGCGCAATTAGCGCAGATTTTGCCACTGCAATTGATTCATCCGGAAGGCTTTGAACTGCTGACCGATGGACCGAAACAGCGTCGGGCTTTTATCGATTGGGGCGTGTTTCATACTGAGCCCGCTTTTTTCGATGCATGGGGGCGGTTTAAGCGGCTCAGCAAGCAGCGTAATGCACTGCTAAAAAGCGCGCAAAGCTATCGTGAACTCAGTTATTGGGATCAAGAATTAGCCCGTTTGGCTGAGCAGATTGACCAGTGGCGTGAAAGTTACGTCAATCAATTGAAAAATGTGGCAGAGCAGTTATGTCGCACATTTTTGCCAGAATTCGATATCGACTTAAAGTATTATCGAGGTTGGGAAAAAGATCAGCCTTATCAATCGATTCTGGAAAAAAACTTCGAACGGGATCAGCAGTTGGGCTATACCTTTAGCGGGCCTAACAAAGCGGATTTGCGGATTAAAGTGAACGCAACCCCGGTGGAAGATGTGTTGTCGCGGGGACAACTGAAGTTGATGGTGTGTGCGTTGCGCGTGGCGCAAGGGCAGCACCTGACCGAGTTGACGGGAAAGCAATGCATTTATCTTATCGATGATTTTGCTTCCGAATTGGATAGCCTACGTCGCCAACGTTTGGCAGATAGCCTAAAAGGGACGGGGGCGCAGGTTTTTGTAAGTTCTATTACCGAAAGCCAAGTGGCGGACATGTTAGATGAGTCCAGTAAAACATTTCACGTTGCCCACGGTGTAATAGAGCAAGGATAA
Filtered Gene (Forward): ATGATCCAACAGTTTCGTAATATTAAAGCCTGTGATATTCGCCTATCAGCAGGCTTTAACTTTCTTATCGGGCCGAATGGCAGCGGTAAAACCAGCGTGCTTGAAGCGATTTATCTGCTCGGCCATGGCCGTTCATTTAAAAGCTCGCTGACCGGACGGATCATTCAAAATGAGTGCTCAGAACTGTTTGTGCATGGCCGGATTTGTGAGCATTCTTTGAGCTCTGATCAATTTGAGCTACCCGTTGGCATTAATAAGCAGCGTGACGGCTCAACTGAGGTTAAAATAGGCGGCCAAACAGGGCAAAAATTGGCGCAATTAGCGCAGATTTTGCCACTGCAATTGATTCATCCGGAAGGCTTTGAACTGCTGACCGATGGACCGAAACAGCGTCGGGCTTTTATCGATTGGGGCGTGTTTCATACTGAGCCCGCTTTTTTCGATGCATGGGGGCGGTTTAAGCGGCTCAGCAAGCAGCGTAATGCACTGCTAAAAAGCGCGCAAAGCTATCGTGAACTCAGTTATTGGGATCAAGAATTAGCCCGTTTGGCTGAGCAGATTGACCAGTGGCGTGAAAGTTACGTCAATCAATTGAAAAATGTGGCAGAGCAGTTATGTCGCACATTTTTGCCAGAATTCGATATCGACTTAAAGTATTATCGAGGTTGGGAAAAAGATCAGCCTTATCAATCGATTCTGGAAAAAAACTTCGAACGGGATCAGCAGTTGGGCTATACCTTTAGCGGGCCTAACAAAGCGGATTTGCGGATTAAAGTGAACGCAACCCCGGTGGAAGATGTGTTGTCGCGGGGACAACTGAAGTTGATGGTGTGTGCGTTGCGCGTGGCGCAAGGGCAGCACCTGACCGAGTTGACGGGAAAGCAATGCATTTATCTTATCGATGATTTTGCTTCCGAATTGGATAGCCTACGTCGCCAACGTTTGGCAGATAGCCTAAAAGGGACGGGGGCGCAGGTTTTTGTAAGTTCTATTACCGAAAGCCAAGTGGCGGACATGTTAGATGAGTCCAGTAAAACATTTCACGTTGCCCACGGTGTAATAGAGCAAGGATAA


## Rosalind 72

ssh Alquwara@ilogin.ibex.kaust.edu.sa

# Create a Directory for the ORF Finder and explore
mkdir ~/orf_finder
cd ~/orf_finder


## Create the Python Script: Use nano to create and edit the orf_finder.py file:

nano orf_finder.py
# the code inside is : 

#!/usr/bin/env python
import sys

# Codon table mapping DNA codons to amino acids
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def flipComp(Text):
    """Compute the reverse complement of a DNA sequence."""
    result = ""
    for letter in Text:
        if letter == 'A':
            result += 'T'
        elif letter == 'T':
            result += 'A'
        elif letter == 'G':
            result += 'C'
        else:
            result += 'G'
    flipped = result[::-1]  # Reverse the sequence
    return flipped

def read_fasta(filename):
    """Reads a single FASTA file and returns the sequence."""
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

def translate_dna(dna_sequence):
    """Translates a DNA sequence into a protein sequence using the codon table."""
    protein = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = CODON_TABLE.get(codon, "")
            if amino_acid == '_':  # Stop codon
                break
            protein += amino_acid
    return protein

def find_orfs(sequence):
    """Finds all ORFs in all six reading frames (both forward and reverse)."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = set()

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = sequence[i:j+3]
                        orfs.add(translate_dna(orf))
                        break

    # Check reverse reading frames
    rev_sequence = flipComp(sequence)
    for frame in range(3):
        for i in range(frame, len(rev_sequence), 3):
            codon = rev_sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(rev_sequence), 3):
                    stop_codon = rev_sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = rev_sequence[i:j+3]
                        orfs.add(translate_dna(orf))
                        break

    return orfs

def main():
    if len(sys.argv) != 2:
        print("Usage: python orf_finder.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

    # Find all ORFs and translate them
    orfs = find_orfs(sequence)

    # Output all distinct protein strings
    for protein in orfs:
        print(protein)

if __name__ == "__main__":
    main()


## Create a Sample FASTA Input File: 
nano input.fasta

# add whatever rosalind give inside 

## Activate Python Environment (if needed):
source activate myenv

# Run the Script: Now, run the ORF finder script with the input FASTA file:
# You should see output printed out from both the forward and reverse reading frames.


python orf_finder.py input.fasta


This is the text file copied from ibex how to use the above functions : 

(myenv) [alquwara@login509-02-l gene_finder]$ cd
(myenv) [alquwara@login509-02-l ~]$ mkdir ~/orf_finder
(myenv) [alquwara@login509-02-l ~]$ cd ~/orf_finder
(myenv) [alquwara@login509-02-l orf_finder]$ nano orf_finder.py
(myenv) [alquwara@login509-02-l orf_finder]$ (myenv) [alquwara@login509-02-l orf_finder]$ nano input.fasta
(myenv) [alquwara@login509-02-l orf_finder]$ (myenv) [alquwara@login509-02-l orf_finder]$ python orf_finder.py
(myenv) [alquwara@login509-02-l orf_finder]$ source activate myenv
(myenv) [alquwara@login509-02-l orf_finder]$ python orf_finder.py
(myenv) [alquwara@login509-02-l orf_finder]$ nano orf_finder.py
(myenv) [alquwara@login509-02-l orf_finder]$ (myenv) [alquwara@login509-02-l orf_finder]$ nano input.fasta
(myenv) [alquwara@login509-02-l orf_finder]$ (myenv) [alquwara@login509-02-l orf_finder]$ source activate myenv
(myenv) [alquwara@login509-02-l orf_finder]$ python orf_finder.py input.fasta
MLLGSFRLIPKETLIQVAGSSPCNLS
MTPRLGLESLLE
M
MGMTPRLGLESLLE

## Task 4: Process All 14 Genomes and Find ORFs


cd
#copy all 14 genome from Reem file to orf_finder

cp ~/gene_finder/Reem/*.fna ~/orf_finder/

#move to the file 
cd ~/orf_finder
#explore 
ls *.fna
 # Output
# GCA_000006745.1_ASM674v1_genomic.fna  GCA_000008545.1_ASM854v1_genomic.fna  GCA_000008745.1_ASM874v1_genomic.fna
# GCA_000006825.1_ASM682v1_genomic.fna  GCA_000008565.1_ASM856v1_genomic.fna  GCA_000008785.1_ASM878v1_genomic.fna
 # GCA_000006865.1_ASM686v1_genomic.fna  GCA_000008605.1_ASM860v1_genomic.fna  GCA_000027305.1_ASM2730v1_genomic.fna
# GCA_000007125.1_ASM712v1_genomic.fna  GCA_000008625.1_ASM862v1_genomic.fna  GCA_000091085.2_ASM9108v2_genomic.fna


## Create a Bash Script
nano run_orf_finder.sh

##  Paste the following code inside nano:
#!/bin/bash

# Loop through all .fna files in the current directory
for genome in *.fna; do
    echo "Processing $genome..."
    # Run the ORF finder and output to a corresponding .txt file
    python ../orf_finder.py $genome > ${genome%.fna}_orfs.txt
    echo "Output saved to ${genome%.fna}_orfs.txt"
done


## Step 2: Make the Script Executable

chmod +x run_orf_finder.sh

#Activate Your Python Environment 
source activate myenv


# Run the Script 
./run_orf_finder.sh

#  Verify the Output
# GCA_000006745.1_ASM674v1_genomic_orfs.txt  GCA_000008605.1_ASM860v1_genomic_orfs.txt
# GCA_000006825.1_ASM682v1_genomic_orfs.txt  GCA_000008625.1_ASM862v1_genomic_orfs.txt
# GCA_000006865.1_ASM686v1_genomic_orfs.txt  GCA_000008725.1_ASM872v1_genomic_orfs.txt
# GCA_000007125.1_ASM712v1_genomic_orfs.txt  GCA_000008745.1_ASM874v1_genomic_orfs.txt
# GCA_000008525.1_ASM852v1_genomic_orfs.txt  GCA_000008785.1_ASM878v1_genomic_orfs.txt
# GCA_000008545.1_ASM854v1_genomic_orfs.txt  GCA_000027305.1_ASM2730v1_genomic_orfs.txt
# GCA_000008565.1_ASM856v1_genomic_orfs.txt  GCA_000091085.2_ASM9108v2_genomic_orfs.txt

#  inspect one of the files:
less GCA_000006745.1_ASM674v1_genomic_orfs.txt 


# Add the output files to Git:

git add *.txt

# Commit the changes:
 
git commit -m "Add ORF results for 14 genomes"


# Push to GitHub:

git push origin master






## Task 5: extend your existing ORF finder to include a filter that discards ORFs shorter than a specified length, with the default being 100 codons.

# Navigate to the correct directory:
cd ~/orf_finder


#Create a new branch for Task 5 (to keep changes separate):
git checkout -b task5-filter-length

## Create a new script for Task 5
nano orf_finder_with_filter.py

#inside this nano this code:

#!/usr/bin/env python
import sys

# Codon table mapping DNA codons to amino acids
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def flipComp(Text):
    """Compute the reverse complement of a DNA sequence."""
    result = ""
    for letter in Text:
        if letter == 'A':
            result += 'T'
        elif letter == 'T':
            result += 'A'
        elif letter == 'G':
            result += 'C'
        else:
            result += 'G'
    flipped = result[::-1]  # Reverse the sequence
    return flipped

def read_fasta(filename):
    """Reads a single FASTA file and returns the sequence."""
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

def translate_dna(dna_sequence):
    """Translates a DNA sequence into a protein sequence using the codon table."""
    protein = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = CODON_TABLE.get(codon, "")
            if amino_acid == '_':  # Stop codon
                break
            protein += amino_acid
    return protein

def find_orfs(sequence, min_length=100):
    """Finds all ORFs in all six reading frames (both forward and reverse) and filters by length."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = set()

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = sequence[i:j+3]
                        if len(orf) >= min_length * 3:  # Filter by length (in nucleotides)
                            orfs.add(translate_dna(orf))
                        break

    # Check reverse reading frames
    rev_sequence = flipComp(sequence)
    for frame in range(3):
        for i in range(frame, len(rev_sequence), 3):
            codon = rev_sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(rev_sequence), 3):
                    stop_codon = rev_sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf = rev_sequence[i:j+3]
                        if len(orf) >= min_length * 3:  # Filter by length (in nucleotides)
                            orfs.add(translate_dna(orf))
                        break

    return orfs

def main():
    if len(sys.argv) != 3:
        print("Usage: python orf_finder_with_filter.py <input_fasta> <min_length>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    min_length = int(sys.argv[2])  # Take the minimum length as input

    sequence = read_fasta(fasta_file)

    # Find all ORFs and filter by length
    orfs = find_orfs(sequence, min_length=min_length)

    output_file = fasta_file.replace('.fna', '_filtered_orfs.txt')

    # Output all distinct protein strings to a file
    with open(output_file, 'w') as f:
        for protein in orfs:
            f.write(protein + "\n")

    print(f"Filtered ORFs saved to {output_file}")

if __name__ == "__main__":
    main()


#Make the script executable:

chmod +x orf_finder_with_filter.py

# Create a directory to store the output for Task 5:

mkdir task5_outputs


# Run the script for all 14 genomes:
nano run_task5_orf_finder.sh


# inside nano this code
#!/bin/bash
for file in *.fna; do
    echo "Processing $file..."
    python orf_finder_with_filter.py "$file" 100
    mv "${file%.fna}_filtered_orfs.txt" task5_outputs/
done
echo "All ORFs have been filtered and saved."


# Make the script executable:

chmod +x run_task5_orf_finder.sh

# Run the script to process all 14 genomes:

./run_task5_orf_finder.sh

#Check the task5_outputs folder for the results:
cd task5_outputs
ls
#output: 
#GCA_000006745.1_ASM674v1_genomic_filtered_orfs.txt  GCA_000008605.1_ASM860v1_genomic_filtered_orfs.txt
#GCA_000006825.1_ASM682v1_genomic_filtered_orfs.txt  GCA_000008625.1_ASM862v1_genomic_filtered_orfs.txt
#GCA_000006865.1_ASM686v1_genomic_filtered_orfs.txt  GCA_000008725.1_ASM872v1_genomic_filtered_orfs.txt
#GCA_000007125.1_ASM712v1_genomic_filtered_orfs.txt  GCA_000008745.1_ASM874v1_genomic_filtered_orfs.txt
#GCA_000008525.1_ASM852v1_genomic_filtered_orfs.txt  GCA_000008785.1_ASM878v1_genomic_filtered_orfs.txt
#GCA_000008545.1_ASM854v1_genomic_filtered_orfs.txt  GCA_000027305.1_ASM2730v1_genomic_filtered_orfs.txt
#GCA_000008565.1_ASM856v1_genomic_filtered_orfs.txt  GCA_000091085.2_ASM9108v2_genomic_filtered_orfs.txt

#Go back to the parent directory:
cd ..

# The script files orf_finder_with_filter.py and run_task5_orf_finder.sh are untracked. include these files in your commit:

git add orf_finder_with_filter.py run_task5_orf_finder.sh

# Add the files from the task5_outputs directory to the Git stage:
git add task5_outputs/*

# commit the changes:
git commit -m "Add filtered ORFs output for Task 5 and Task 5 scripts"

# push files: 
git push origin task5-filter-length







## Task 6 :  implement a ribosome binding site (RBS) filter to identify ORFs with a Shine-Dalgarno sequence (AGGAGG) located 4-20bp upstream of the start codon. This process involves:

 ## Set up the directory for Task 6
cd ~/orf_finder/
mkdir task6_filter_rbs_outputs

# Copy the Genome Files
#Go back to the parent directory:
cd ~/gene_finder/Reem/

cp *.fna ~/orf_finder/task6_filter_rbs_outputs/

# Copy the previous script for modification:
cp orf_finder_with_filter.py orf_finder_with_rbs.py

# Create the Python Script for Task 6

cd ~/orf_finder/

#Copy the existing ORF finder script and name it for Task 6:
cp orf_finder_with_filter.py orf_finder_with_rbs_filter.py

#Open the new script 
nano orf_finder_with_rbs_filter.py

#inside nano : 
#!/usr/bin/env python
import sys

# Codon table mapping DNA codons to amino acids
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

# RBS sequence and upstream distance for Shine-Dalgarno sequence
RBS_SEQUENCE = "AGGAGG"
UPSTREAM_DISTANCE = 20

# Function to find reverse complement
def flipComp(sequence):
    """Compute the reverse complement of a DNA sequence, handling ambiguous nucleotide codes."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N'
    }
    return ''.join([complement.get(base, 'N') for base in sequence[::-1]])

# Function to read FASTA file
def read_fasta(filename):
    """Reads a single FASTA file and returns the sequence."""
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines[1:]])  # Ignore the first line (FASTA header)
    return sequence

# Function to translate DNA into a protein sequence
def translate_dna(dna_sequence):
    """Translates a DNA sequence into a protein sequence."""
    protein = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = CODON_TABLE.get(codon, "")
            if amino_acid == '_':  # Stop codon
                break
            protein += amino_acid
    return protein

# Function to check for RBS sequence
def find_rbs(sequence, start_index):
    """Check for RBS sequence within a specified upstream distance of the start codon."""
    upstream_region = sequence[max(0, start_index - UPSTREAM_DISTANCE):start_index]
    return RBS_SEQUENCE in upstream_region

# Function to find ORFs
def find_orfs(sequence, min_length=100):
    """Finds all ORFs in all six reading frames (both forward and reverse)."""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []

    # Check forward reading frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(sequence), 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        if find_rbs(sequence, i) and (j - i + 3) // 3 >= min_length:
                            orf = sequence[i:j+3]
                            orfs.append(translate_dna(orf))
                        break

    # Check reverse reading frames
    rev_sequence = flipComp(sequence)
    for frame in range(3):
        for i in range(frame, len(rev_sequence), 3):
            codon = rev_sequence[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(rev_sequence), 3):
                    stop_codon = rev_sequence[j:j+3]
                    if stop_codon in stop_codons:
                        if find_rbs(rev_sequence, i) and (j - i + 3) // 3 >= min_length:
                            orf = rev_sequence[i:j+3]
                            orfs.append(translate_dna(orf))
                        break

    return orfs

# Main function to execute the ORF finder
def main():
    if len(sys.argv) != 2:
        print("Usage: python orf_finder_with_rbs_filter.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)

    # Find ORFs and filter them
    orfs = find_orfs(sequence)

    # Write the ORFs to output file
    output_file = f"./task6_filter_rbs_outputs/{fasta_file.split('/')[-1].replace('.fna', '_task6_filtered_orfs.txt')}"
    with open(output_file, 'w') as f:
        for orf in orfs:
            f.write(f"{orf}\n")

    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()


##  Create the run_task6_orf_finder.sh script:
nano run_task6_orf_finder.sh

#inside nano this code: 
#!/bin/bash

# Make sure the output directory exists
mkdir -p task6_filter_rbs_outputs

# Loop through all FASTA files and run the ORF finder script
for fasta_file in ../Reem/*.fna; do
    echo "Processing $fasta_file..."
    python orf_finder_with_rbs_filter.py "$fasta_file"
done

## Prepare and Run the Scripts

cd ~/orf_finder


#  Make the shell script executable

chmod +x run_task6_orf_finder.sh

mkdir -p task6_filter_rbs_outputs

# Run the script
./run_task6_orf_finder.sh


# Check the results:
ls task6_filter_rbs_outputs


# Add the results to Git:

git add task6_filter_rbs_outputs/*


# Commit the changes:

git commit -m "Add filtered ORF output for Task 6 with RBS filtering"



# Push the changes to a new branch:
git push origin master

