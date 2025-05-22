#!/usr/bin/env python3
import sys
import re

def process_line(line, ids):
    fields = line.strip().split("\t")
    sum_value = sum(float(field) for field in fields[1:])
    return f"{ids.pop(0)}\t{line.strip()}\t{sum_value}"

def main():
    # Check if the number of arguments is correct
    if len(sys.argv) != 4:
        print("Usage: {} fasta_file table_file output_file".format(sys.argv[0]))
        sys.exit(1)

    # Read IDs from the FASTA file
    with open(sys.argv[1], 'r') as fasta_file:
        fasta_ids = re.findall(r'^>(\w+)', fasta_file.read(), re.MULTILINE)

    # Read the table and add the IDs as the first column
    with open(sys.argv[2], 'r') as table_file:
        table_lines = table_file.readlines()

    # Write the output file with IDs and sum
    with open(sys.argv[3], 'w') as output_file:
        output_file.write("id\t{}\tsum\n".format(table_lines[0].strip()))
        for line, fasta_id in zip(table_lines[1:], fasta_ids):
            output_file.write(process_line(line, [fasta_id]) + '\n')

if __name__ == "__main__":
    main()
