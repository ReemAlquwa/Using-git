

## Task 1 week 5
# create a directory specifically for Week 5
cd ~/gene_finder
mkdir week5_genome_annotation
cd week5_genome_annotation

# README for Task Solutions
nano README.md
nano task1_amino_acid_count.py


#inside nano add : 

# Define the amino acid sequence
amino_acid_sequence = "KVRMFTSELDIMLSVNGPADQIKYFCRHWT"

# Count the number of amino acids (excluding stop codon)
num_amino_acids = len(amino_acid_sequence)

# Calculate the total number of bases (including stop codon)
total_bases = (num_amino_acids * 3) + 3

# Print the results
print(f"Number of amino acids: {num_amino_acids}")
print(f"Total number of bases (including stop codon): {total_bases}")



# Count Amino Acids and Bases in ORF

 (myenv) [alquwara@login509-02-l ~]$ python task1_amino_acid_count.py

 ##output:
# Number of amino acids: 30
# Total number of bases (including stop codon): 93



cd ~/gene_finder/week5_genome_annotation

git status

git add task1_amino_acid_count.py README.md 'AA&BASES_README.md' 

git commit -m "Add Week 5 task scripts and README files"

git push origin master




## Task 2 week 5

#Find the genome file:

find /home/alquwara -name "GCA_000006745.1_ASM674v1_genomic.fna"


#output: 
[alquwara@login509-02-l ~]$ find /home/alquwara -name "GCA_000006745.1_ASM674v1_genomic.fna"
/home/alquwara/ncbi_dataset/data/Reem/GCA_000006745.1_ASM674v1_genomic.fna

#Copy the genome file to your Week 5 directory:

cp /home/alquwara/ncbi_dataset/data/Reem/GCA_000006745.1_ASM674v1_genomic.fna ~/gene_finder/week5_genome_annotation/

# After copying, navigate to the week5_genome_annotation directory:

cd ~/gene_finder/week5_genome_annotation

##Run Prodigal on this genome file to annotate genes:

#Use module command to check for Prodigal:
module avail prodigal

#output: 
[alquwara@login509-02-l week5_genome_annotation]$ module avail prodigal
---------------------------------------------------------------------------- /sw/rl9c/modulefiles/applications -----------------------------------------------------------------------------
prodigal/2.6.3

Key:
modulepath


#load it : 

module load prodigal/2.6.3

# run the Prodigal command on the genome:

prodigal -i GCA_000006745.1_ASM674v1_genomic.fna -o GCA_000006745.1_ASM674v1_genomic_output.txt -a GCA_000006745.1_ASM674v1_genomic_proteins.faa

#output: 
[alquwara@login509-02-l week5_genome_annotation]$ prodigal -i GCA_000006745.1_ASM674v1_genomic.fna -o GCA_000006745.1_ASM674v1_genomic_output.txt -a GCA_000006745.1_ASM674v1_genomic_proteins.faa
-------------------------------------
PRODIGAL v2.6.3 [February, 2016]
Univ of Tenn / Oak Ridge National Lab
Doug Hyatt, Loren Hauser, et al.
-------------------------------------
Request:  Single Genome, Phase:  Training
Reading in the sequence(s) to train...4033488 bp seq created, 47.49 pct GC
Locating all potential starts and stops...219838 nodes
Looking for GC bias in different frames...frame bias scores: 2.21 0.16 0.62
Building initial set of genes to train from...done!
Creating coding model and scoring nodes...done!
Examining upstream regions and training starts...done!
-------------------------------------
Request:  Single Genome, Phase:  Gene Finding
Finding genes in sequence #1 (2961149 bp)...done!
Finding genes in sequence #2 (1072315 bp)...done!

#This will generate the output files:
# .txt for the annotated genes
and 
# .faa for the protein sequences

# Count the gene entries in the .txt output file:
 grep -c "CDS" GCA_000006745.1_ASM674v1_genomic_output.txt

 ## output
## 3594

#or 
#Count the number of sequences in the protein file:

grep -c "^>" GCA_000006745.1_ASM674v1_genomic_proteins.faa

##output:
## 3594




## Task 3 week 5

#Go to directory:
 cd /home/alquwara/ncbi_dataset/data/Reem

#list  genomes:
 ls

 #output:
GCA_000006745.1_ASM674v1_genomic.fna  GCA_000008525.1_ASM852v1_genomic.fna  GCA_000008625.1_ASM862v1_genomic.fna  GCA_000027305.1_ASM2730v1_genomic.fna
GCA_000006825.1_ASM682v1_genomic.fna  GCA_000008545.1_ASM854v1_genomic.fna  GCA_000008725.1_ASM872v1_genomic.fna  GCA_000091085.2_ASM9108v2_genomic.fna
GCA_000006865.1_ASM686v1_genomic.fna  GCA_000008565.1_ASM856v1_genomic.fna  GCA_000008745.1_ASM874v1_genomic.fna
GCA_000007125.1_ASM712v1_genomic.fna  GCA_000008605.1_ASM860v1_genomic.fna  GCA_000008785.1_ASM878v1_genomic.fna

##Create a Shell Script:
# go to the directory week5_genome_annotation:
cd ~/gene_finder/week5_genome_annotation

#Create a shell script file:
nano run_prodigal_all_genomes.sh

#Inside this shell script, write the following code:

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


#Make the script executable:
chmod +x run_prodigal_all_genomes.sh

#Run the script:
./run_prodigal_all_genomes.sh

#Check the Output:
#View the gene_counts.txt file to verify the counts:
less gene_counts.txt

#Commit and Push the Results to GitHub:
git add run_prodigal_all_genomes.sh gene_counts.txt

git commit -m "Add Task 3 solution for Week 5 - Prodigal run on all genomes"

git push origin master


## Task 4 week 5

#working directory:
cd ~/gene_finder/week5_genome_annotation

#Create a new shell script using nano:

nano run_prokka_all_genomes.sh

#inside nano: 
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

#Run the updated script:
chmod +x run_prokka_all_genomes.sh
./run_prokka_all_genomes.sh


#Check the prokka_gene_counts.txt file:
less prokka_gene_counts.txt

#output:

Genome File - Annotated CDS Count by Prokka
GCA_000006745.1_ASM674v1_genomic.fna - 3589
GCA_000006825.1_ASM682v1_genomic.fna - 2028
GCA_000006865.1_ASM686v1_genomic.fna - 2383
GCA_000007125.1_ASM712v1_genomic.fna - 3150
GCA_000008525.1_ASM852v1_genomic.fna - 1577
GCA_000008545.1_ASM854v1_genomic.fna - 1861
GCA_000008565.1_ASM856v1_genomic.fna - 3245
GCA_000008605.1_ASM860v1_genomic.fna - 1001
GCA_000008625.1_ASM862v1_genomic.fna - 1771
GCA_000008725.1_ASM872v1_genomic.fna - 892
GCA_000008745.1_ASM874v1_genomic.fna - 1058
GCA_000008785.1_ASM878v1_genomic.fna - 1504
GCA_000027305.1_ASM2730v1_genomic.fna - 1748
GCA_000091085.2_ASM9108v2_genomic.fna - 1056



#Navigate to the project directory :
cd ~/gene_finder/week5_genome_annotation

#Check the status to see which files are untracked or modified::

git status

#Push files to GitHub

git add prokka_gene_counts.txt run_prokka_all_genomes.sh prokka_output/ GCA_*.txt

git commit -m "Add Week 5 Task 4: Prokka genome annotation results and shell script"

git push origin master



## Task 5 week 5


#go to directory
[alquwara@login509-02-l week5_genome_annotation]$ cd ~/gene_finder/week5_genome_annotation
t_unique_genes.sh

#open shell script (nano)
[alquwara@login509-02-l week5_genome_annotation]$ nano extract_unique_genes.sh

#Make the script executable
[alquwara@login509-02-l week5_genome_annotation]$ chmod +x extract_unique_genes.sh


[alquwara@login509-02-l week5_genome_annotation]$ ./extract_unique_genes.sh

#Run the Script
#Now, run the script to extract unique gene names from all the .gff files in the prokka_output directory:

First five unique gene names:
aaaT
aaeA
aaeA_1
aaeA_2
aaeB

#Commit the Results and Script to Git
git add extract_unique_genes.sh unique_gene_names.txt

git commit -m "Add script for extracting unique gene names annotated by Prokka"

git push origin master



