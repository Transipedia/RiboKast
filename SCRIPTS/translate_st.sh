#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 4 ]; then
    echo "Usage: $0 input_file output_file temp_fasta temp_result"
    exit 1
fi

# Assign arguments to variables
input_file="$1"
output_file="$2"
temp_fasta="$3"
temp_result="$4"

# Read the first line of the input file (the header)
IFS=$'\t' read -r header_line < "$input_file"

# Write the header to the output file, adding the "translated seq" column
echo -e "$header_line\ttranslated_seq" > "$output_file"

# Read the rest of the file line by line, ignoring the header
tail -n +2 "$input_file" | while IFS=$'\t' read -r seq ID_contig P0 P1 P2 Dominant_Phase phase_plus_1; do
    # Write the sequence to a temporary FASTA file
    echo ">seq" > "$temp_fasta"
    echo "$seq" >> "$temp_fasta"
    
    # Use the value of "Phase+1" directly as the translation phase
    phase="$phase_plus_1"
    
    # Translate the sequence into protein with seqkit according to the phase specified in "Phase+1"
    seqkit translate -F -f "$phase" "$temp_fasta" --line-width 7000 > "$temp_result"
    
    # Retrieve the translated protein sequence (all lines except the first one)
    translated_sequence=$(sed -n '2p' "$temp_result")
    echo "$translated_sequence"

    
    # Write the translated sequence to the output file
    echo -e "$seq\t$ID_contig\t$P0\t$P1\t$P2\t$Dominant_Phase\t$phase_plus_1\t$translated_sequence" >> "$output_file"
    
done

echo "The sequences have been successfully translated. See the file: $output_file"
