
# Week3 All tasks



 #input
 scp\Users\Alquwara\Downloads\ncbi.zip alquwara@ilogin.ibex.kaust.edu.sa:/home/alquwara

#output

 ncbi.zip                                                                              100% 8192KB   7.1MB/s   00:01

#input

 unzip ncbi.zip

 #output

 Archive:  ncbi.zip
  inflating: README.md
  inflating: ncbi_dataset/data/data_summary.tsv
  inflating: ncbi_dataset/data/assembly_data_report.jsonl
  inflating: ncbi_dataset/data/GCA_000006745.1/GCA_000006745.1_ASM674v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000006825.1/GCA_000006825.1_ASM682v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000006865.1/GCA_000006865.1_ASM686v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000007125.1/GCA_000007125.1_ASM712v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008545.1/GCA_000008545.1_ASM854v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008565.1/GCA_000008565.1_ASM856v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008605.1/GCA_000008605.1_ASM860v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008625.1/GCA_000008625.1_ASM862v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008725.1/GCA_000008725.1_ASM872v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008785.1/GCA_000008785.1_ASM878v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008525.1/GCA_000008525.1_ASM852v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008745.1/GCA_000008745.1_ASM874v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000027305.1/GCA_000027305.1_ASM2730v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000091085.2/GCA_000091085.2_ASM9108v2_genomic.fna
  inflating: ncbi_dataset/data/dataset_catalog.json
  inflating: md5sum.txt

  #input
  



## Copy file and unzip


 #input
 scp\Users\Alquwara\Downloads\ncbi.zip alquwara@ilogin.ibex.kaust.edu.sa:/home/alquwara

#output

 ncbi.zip                                                                              100% 8192KB   7.1MB/s   00:01

#input

 unzip ncbi.zip

 #output

 Archive:  ncbi.zip
  inflating: README.md
  inflating: ncbi_dataset/data/data_summary.tsv
  inflating: ncbi_dataset/data/assembly_data_report.jsonl
  inflating: ncbi_dataset/data/GCA_000006745.1/GCA_000006745.1_ASM674v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000006825.1/GCA_000006825.1_ASM682v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000006865.1/GCA_000006865.1_ASM686v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000007125.1/GCA_000007125.1_ASM712v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008545.1/GCA_000008545.1_ASM854v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008565.1/GCA_000008565.1_ASM856v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008605.1/GCA_000008605.1_ASM860v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008625.1/GCA_000008625.1_ASM862v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008725.1/GCA_000008725.1_ASM872v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008785.1/GCA_000008785.1_ASM878v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008525.1/GCA_000008525.1_ASM852v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000008745.1/GCA_000008745.1_ASM874v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000027305.1/GCA_000027305.1_ASM2730v1_genomic.fna
  inflating: ncbi_dataset/data/GCA_000091085.2/GCA_000091085.2_ASM9108v2_genomic.fna
  inflating: ncbi_dataset/data/dataset_catalog.json
  inflating: md5sum.txt
## smallest genome in data

#input
tail -n +2 data_summary.tsv | sort -t$'\t' -k11 -n | head -n 1

#output

Chlamydia trachomatis D/UW-3/CX         strain: D/UW-3/CX       272561  ASM872v1        GCA_000008725.1 GenBank Annotation submitted by ChGP    Complete Genome 1042519 1042519 2001-01-09      939     PRJNA45 SAMN02603114
## largest

#input
sort -t$'\t' -k11 -n data_summary.tsv | tail -n 1


#output
Vibrio cholerae O1 biovar El Tor str. N16961            strain: N16961  243277  ASM674v1        GCA_000006745.1 GenBank                                                                     Annotation submitted by TIGR     Complete Genome 2961149 4033464 2001-01-09      4007    PRJNA36 SAMN02603969

## genome size smallest and largest

#input (small genome size)

 tail -n +2 data_summary.tsv | sort -t$'\t' -k11 -n | head -n 1 | cut -f11

#output 

1042519

#input (large genome size)
sort -t$'\t' -k11 -n data_summary.tsv | tail -n 1 | cut -f11


#output

4033464
## Find the number of genomes that contain at least two “c” in the species name. How many of the species names contain two or more “c” but do not contain the word “coccus”?

#input for number of genomes that contain at least two “c” in the species name:
grep -E "\b\w*c\w*c\w*\b" *.fna | wc -l


#output  number of genomes that contain at least two “c” in the species name:
7

#input How many of the species names contain two or more “c” but do not contain the word “coccus”? 

 grep -E "\b\w*c\w*c\w*\b" *.fna | grep -v "coccus" | wc -l
#output How many of the species names contain two or more “c” but do not contain the word “coccus”?
2




## Use the find command to find all genome files (FASTA) larger than 3 megabyte. How many are there


#input
find . -name "*.fna" -size +3M | wc -l

#output
3