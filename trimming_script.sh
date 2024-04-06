#!/bin/bash

# Path to Trimmomatic JAR file
TRIMMOMATIC="trimmomatic.jar"

# Loop over pairs of input files
for ((i=1; i<=$#; i+=2)); do
    input_forward="${!i}"
    input_reverse="${!((i+1))}"
    
    # Output file names
    output_forward_paired="trimmed_${input_forward}"
    output_forward_unpaired="unpaired_${input_forward}"
    output_reverse_paired="trimmed_${input_reverse}"
    output_reverse_unpaired="unpaired_${input_reverse}"

    # Run Trimmomatic
    java -jar "$TRIMMOMATIC" PE -phred33 "$input_forward" "$input_reverse" \
        "$output_forward_paired" "$output_forward_unpaired" \
        "$output_reverse_paired" "$output_reverse_unpaired" \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
