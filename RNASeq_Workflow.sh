# Download and install SRA toolkit
# follow the installation instructions from the official SRA toolkit website

# Download SRA accessions
prefetch SRR15838464 SRR15838465 SRR15838468 SRR15838469

# Convert downloaded SRA files to FASTQ format
fastq-dump --gzip --split-files SRR15838464 SRR15838465 SRR15838468 SRR15838469

# Create a directory for FastQC results
mkdir /mnt/f/sra/fastqc_results

# Run FastQC on the downloaded FASTQ files and save results to the directory
fastqc /mnt/f/sra/*.fastq.gz -o /mnt/f/sra/fastqc_results

# Write the shell script for trimming (only trimming is needed)
# Run the shell script for trimming
./trimming_script.sh SRR15838464_1.fastq SRR15838464_2.fastq SRR15838465_1.fastq SRR15838465_2.fastq SRR15838468_1.fastq SRR15838468_2.fastq SRR15838469_1.fastq SRR15838469_2.fastq

# Build an index for the reference genome
hisat2-build IRGSP-1.0_genome.fasta IRGSP-1.0_index

# Unzip the downloaded FASTQ files
gunzip *.fastq.gz

# Align paired-end reads to the reference genome using HISAT2
hisat2 -x IRGSP-1.0_index -1 SRR15838464_1.fastq -2 SRR15838464_2.fastq -S aligned_reads_SRR15838464.sam
# Repeat this command for each pair of FASTQ files, replacing the filenames and output SAM file names accordingly

# Convert SAM to BAM and sort
samtools view -bS aligned_reads_SRR15838464.sam | samtools sort -o aligned_reads_sorted_SRR15838464.bam
# Repeat this command for each SAM file, replacing the filenames and output BAM file names accordingly

# Index the sorted BAM files
samtools index aligned_reads_sorted_SRR15838464.bam
# Repeat this command for each sorted BAM file

# Quantify gene expression using FeatureCounts
featureCounts -p -a IRGSP-1.0_representative_transcript_exon_2024-01-11.gtf -o counts_data.txt -T 4 aligned_reads_sorted_SRR15838464.bam aligned_reads_sorted_SRR15838465.bam aligned_reads_sorted_SRR15838468.bam aligned_reads_sorted_SRR15838469.bam
