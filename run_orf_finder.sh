#!/bin/bash

# Loop through all .fna files in the current directory
for genome in *.fna; do
    echo "Processing $genome..."
    # Run the ORF finder and output to a corresponding .txt file
    python orf_finder.py $genome > ${genome%.fna}_orfs.txt
    echo "Output saved to ${genome%.fna}_orfs.txt"
done

