#!/bin/bash

# Load Prodigal module
module load prodigal/2.6.3

# Path to genome files
GENOME_DIR="/home/alquwara/ncbi_dataset/data/Reem"

# Create an output file to store results
echo "Genome File - Annotated Genes Count" > gene_counts.txt

# Loop through all .fna files and run Prodigal
for genome in $GENOME_DIR/*.fna; do
  # Get the genome file name without path
  genome_name=$(basename $genome)

  # Run Prodigal and output to individual files
  prodigal -i $genome -o ${genome_name}_output.txt -a ${genome_name}_proteins.faa

  # Count the number of CDS annotations (genes)
  gene_count=$(grep -c "CDS" ${genome_name}_output.txt)

  # Save the result in gene_counts.txt
  echo "$genome_name - $gene_count" >> gene_counts.txt

  echo "Processed $genome_name with $gene_count genes"
done

# Find the genome with the highest number of genes
echo "Genome with the highest number of genes:"
sort -k3 -n gene_counts.txt | tail -n 1

