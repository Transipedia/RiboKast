import sys

# Function to generate k-mers from a given sequence
def generate_kmers(sequence, k):
    kmers = []
    for i in range(0, len(sequence) - k + 1, 1):
        kmers.append(sequence[i:i+k])
    return kmers

# Function to write k-mers to an output file
def write_kmers_to_file(sequences, output_file, k):
    try:
        with open(output_file, 'w') as f:
            for sequence_id, sequence in sequences.items():
                kmers = generate_kmers(sequence, k)
                for i, kmer in enumerate(kmers, start=1):
                    f.write(f">{sequence_id}_kmer_{i}\n{kmer}\n")
        print("K-mers have been written to the file:", output_file)
    except Exception as e:
        print("An error occurred while writing to the file:", str(e))

# Main script
if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py input_file output_file k")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if the value of k is an integer
    try:
        k = int(sys.argv[3])
    except ValueError:
        print("The value of k must be an integer.")
        sys.exit(1)
    
    try:
        sequences = {}
        # Read the input file and parse sequences
        with open(input_file, 'r') as f:
            current_sequence_id = None
            current_sequence = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_sequence_id is not None:
                        sequences[current_sequence_id] = current_sequence
                    current_sequence_id = line[1:]  # Remove '>' to get the sequence ID
                    current_sequence = ""
                else:
                    current_sequence += line
            if current_sequence_id is not None:
                sequences[current_sequence_id] = current_sequence
        
        # Write k-mers to the output file
        write_kmers_to_file(sequences, output_file, k)
    except FileNotFoundError:
        print("The specified input file was not found.")
    except Exception as e:
        print("An error occurred:", str(e))
