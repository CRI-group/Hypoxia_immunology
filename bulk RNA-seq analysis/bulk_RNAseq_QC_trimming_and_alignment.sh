#!/bin/bash -l
#SBATCH -J analyse_bulkRNAseq
#SBATCH --output=analyse_bulkRNAseq.out
#SBATCH --error=analyse_bulkRNAseq.err
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem= 

#Script to run QC, trim and align bulk RNA seq samples. Made by Iina Koivisto in 2024.
#Make "STAR" and "fastqc" named folders in your working directory prior to starting the analysis.

#Change file names here to your file names
file1="M_19_74_24h_EKRN230069226-1A_22H3C7LT3_L2_1"
file2="M_19_74_24h_EKRN230069226-1A_22H3C7LT3_L2_2"
prefix="M_19_74_24h_"

echo "Starting the analysis of $prefix!"

#Quality control
module load compbio/fastqc/0.11.9
fastqc -o fastqc ${file1}.fq.gz ${file2}.fq.gz
echo "Fastqc done, starting trimming!"
 
#Trimming adapter sequences and reads with quality below 20
module load compbio/trimgalore/0.6.7
trim_galore --paired ${file1}.fq.gz ${file2}.fq.gz
echo "Trimming done, starting fastqc!"

#second QC to see if trimming succeeded
fastqc -o fastqc ${file1}_val_1.fq.gz ${file2}_val_2.fq.gz
echo "Fastqc done, starting alignment!"

#Alignment with STAR compbio/STAR/2.7.10b index created by Anja Hartewig. Located in /lustre/compbio/pub/references/Gencode_v43/STAR_index
module load compbio/STAR/2.7.10b
STAR --genomeDir /lustre/compbio/pub/references/Gencode_v43/STAR_index/ --runThreadN 2 --readFilesIn ${file1}_val_1.fq.gz ${file2}_val_2.fq.gz --outFileNamePrefix STAR/${prefix} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
echo "Alignment done, starting samtools to make index!"

#making bai index using samtools compbio/samtools/1.15.1
module load compbio/samtools/1.15.1
samtools index STAR/${prefix}Aligned.sortedByCoord.out.bam

#echo "Index done, starting fastqc!"

#qc to see how alignment succeeded
#fastqc -o fastqc STAR/${prefix}Aligned.sortedByCoord.out.bam

echo "Analysis ready!"
