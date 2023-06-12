#!/bin/bash

input_file="$1"
output_folder="$2"
output_prefix="output"
seq_count=0
file_count=1
sequence=""
in_sequence=false

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Function to handle opening input file (with or without decompression)
open_input_file() {
    if [[ "$input_file" == *.gz ]]; then
        zcat "$input_file"
    else
        cat "$input_file"
    fi
}

# Iterate through the input file and split into multiple output files
open_input_file | while IFS= read -r line; do
    if [[ $line =~ ^\> ]]; then
        if [ "$in_sequence" = true ]; then
            seq_count=$((seq_count + 1))
            sequence+=$'\n'

            # Check if we need to start a new output file
            if [ $seq_count -eq 1 ]; then
                output_file="${output_folder}/${output_prefix}_${file_count}.fasta"
                echo -e "$sequence" > "$output_file"
            else
                echo -e "$sequence" >> "$output_file"
            fi

            # Check if we reached the millionth sequence
            if [ $seq_count -eq 1000000 ]; then
                seq_count=0
                file_count=$((file_count + 1))
            fi

            sequence=""
        fi

        sequence+="$line"
        in_sequence=true
    else
        sequence+="$line"
    fi
done

# Write the last sequence to the output file
if [ "$in_sequence" = true ]; then
    seq_count=$((seq_count + 1))
    sequence+=$'\n'
    echo -e "$sequence" >> "$output_file"
fi

echo "Split $input_file into $file_count output files."

