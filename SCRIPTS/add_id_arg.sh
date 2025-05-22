#!/bin/bash

# Check the number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <file1> <file2> <output_file>"
    exit 1
fi

# Input and output files from arguments
file1="$1"
file2="$2"
output_file="$3"

# Create an empty output file
> $output_file

# Read file 2 and store peptide IDs and corresponding sequences in an associative array
declare -A peptide_map

while read -r line; do
    if [[ $line == ">"* ]]; then
        current_peptide_id="${line:1}"  # Remove the '>' character
    else
        peptide_map["$line"]="$current_peptide_id"
    fi
done < $file2

# Process file 1 line by line
while IFS=$'\t' read -r sequence rest; do
    # Look for a match in the associative array
    peptide_id="${peptide_map["$sequence"]}"
    
    # If a match is found, add the ID to the line, otherwise put "unknown"
    if [[ -n $peptide_id ]]; then
        echo -e "$peptide_id\t$sequence\t$rest" >> $output_file
    else
        echo -e "ID\t$sequence\t$rest" >> $output_file
    fi
done < $file1

echo "The file $output_file has been created successfully."
