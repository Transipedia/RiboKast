import argparse

def reverse_complement(sequence):
    """Generate the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.translate(complement)[::-1]

def process_fasta(file_path):
    """Process the FASTA file to add reverse complement for _uns IDs."""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    output_lines = []
    for i in range(0, len(lines), 2):
        header = lines[i].strip()
        sequence = lines[i + 1].strip()

        output_lines.append(header + "Seq")
        output_lines.append(sequence)

        if "_uns" in header:
            revcomp_sequence = reverse_complement(sequence)
            output_lines.append(header + "Rev")
            output_lines.append(revcomp_sequence)

    return output_lines

def write_output(output_path, output_lines):
    """Write the processed lines to a new file."""
    with open(output_path, 'w') as f:
        f.write("\n".join(output_lines) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTA file and add reverse complement for _uns IDs.")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_file", help="Path to the output FASTA file")
    args = parser.parse_args()

    # Process the input file and write the output
    output_lines = process_fasta(args.input_file)
    write_output(args.output_file, output_lines)

    print(f"Processed output written to {args.output_file}")
