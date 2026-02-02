import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Add a contig column to the table file.")
parser.add_argument("fasta_file", help="Path to the FASTA file containing contig ID <-> contig mappings")
parser.add_argument("table_file", help="Path to the table file")
parser.add_argument("output_file", help="Path to the output file")

# Parse the arguments
args = parser.parse_args()

# Read the FASTA file and store the contig ID <-> contig mappings
fasta_contig_mapping = {}
with open(args.fasta_file, "r") as fasta_file:
    contig_id = None
    for line in fasta_file:
        if line.startswith(">"):
            contig_id = line.strip()[1:]  # Remove ">" from the contig ID
        else:
            fasta_contig_mapping[contig_id] = line.strip()

# Add the contig column to the table file
with open(args.table_file, "r") as table_file:
    with open(args.output_file, "w") as output_file:
        header = table_file.readline().strip()
        output_file.write("contig\t" + header + "\n")  # Add a new column for the contig identifiers
        for line in table_file:
            contig_id = line.split()[0]  # Extract the contig ID from the table file
            #contig_id = line.split("_", 2)[0]

            print(contig_id)
            if contig_id in fasta_contig_mapping:  # Check if the contig ID exists in the FASTA file
                corresponding_contig = fasta_contig_mapping[contig_id]
                output_file.write(corresponding_contig + "\t" + line)  # Add the corresponding contig ID as the first column
