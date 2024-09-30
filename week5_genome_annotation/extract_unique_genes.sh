#!/bin/bash

# Output file to store unique gene names
echo "Unique Gene Names from Prokka Outputs" > unique_gene_names.txt

# Loop through all Prokka .gff files and extract gene names
for file in prokka_output/*/*.gff; do
  # Extract gene names and append them to the unique_gene_names.txt file
  grep "gene=" $file | sed 's/.*gene=\([^;]*\);.*/\1/' >> all_genes.txt
done

# Sort and remove duplicates to get unique gene names
sort all_genes.txt | uniq > unique_gene_names.txt

# Show the first five unique gene names
echo "First five unique gene names:"
head -n 5 unique_gene_names.txt

