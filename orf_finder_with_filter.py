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

