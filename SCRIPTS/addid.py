import pandas as pd
import sys

# Get command-line arguments
fasta_file = sys.argv[1]
data_file = sys.argv[2]
output_file = sys.argv[3]

# Extract IDs from the FASTA file
ids = []
with open(fasta_file, 'r') as fasta:
    for line in fasta:
        if line.startswith('>'):
            ids.append(line.strip()[1:])

# Read the data file into a DataFrame
data = pd.read_csv(data_file, sep='\t')

# Ensure the number of IDs matches the number of rows in the data
if len(ids) == len(data):
    data.insert(0, 'ID', ids)  # Insert the IDs as the first column
else:
    print("Error: The number of IDs does not match the number of rows in the data file.")
    sys.exit(1)

# Save the modified DataFrame to the output file
data.to_csv(output_file, sep='\t', index=False)
