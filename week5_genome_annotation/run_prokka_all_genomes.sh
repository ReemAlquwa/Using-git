#!/bin/bash

# Path to genome files
GENOME_DIR="/home/alquwara/ncbi_dataset/data/Reem"
PROKKA_OUTPUT_DIR="./prokka_output"

# Create an output file to store results
echo "Genome File - Annotated CDS Count by Prokka" > prokka_gene_counts.txt

# Loop through all .fna files and run Prokka
for genome in $GENOME_DIR/*.fna; do
  # Get the genome file name without path
  genome_name=$(basename $genome .fna)

  # Run Prokka and output to the prokka_output directory
  prokka --force --outdir $PROKKA_OUTPUT_DIR/$genome_name --prefix $genome_name $genome

  # Extract the CDS count from the Prokka output summary
  CDS_count=$(grep "CDS:" $PROKKA_OUTPUT_DIR/$genome_name/$genome_name.txt | awk '{print $2}')

  # Save the result in prokka_gene_counts.txt
  echo "$genome_name.fna - $CDS_count" >> prokka_gene_counts.txt

  echo "Processed $genome_name with $CDS_count CDS"
done

# Find the genome with the highest number of CDS
echo "Genome with the highest number of CDS annotations:"
sort -k3 -n prokka_gene_counts.txt | tail -n 1

