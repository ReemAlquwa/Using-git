#!/bin/bash
for file in *.fna; do
    echo "Processing $file..."
    python orf_finder_with_filter.py "$file" 100
    mv "${file%.fna}_filtered_orfs.txt" task5_outputs/
done
echo "All ORFs have been filtered and saved."

