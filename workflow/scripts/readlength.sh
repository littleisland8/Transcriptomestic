#!/bin/bash

# Provide the path to your FASTQ file
fastq_file="your_fastq_file.fastq"

# Extract sequence lengths and compute total length
total_length=$(awk 'NR%4==2{sum+=length($0)}END{print sum}' "$fastq_file")

# Count the number of reads
read_count=$(grep -c "^@" "$fastq_file")

# Calculate average read length
average_length=$(echo "scale=2; $total_length / $read_count" | bc)

echo "Total bases: $total_length"
echo "Total reads: $read_count"
echo "Average read length: $average_length"
