import sys

def select_peptide_before_first_stop(sequence):
    """
    Extracts the peptide sequence before the first stop codon (*).

    Args:
        sequence (str): Protein sequence.

    Returns:
        tuple: A tuple containing (peptide, start, end) where:
            - peptide (str): The peptide sequence before the first stop codon.
            - start (int): The starting position of the peptide in the sequence (1-indexed).
            - end (int): The ending position of the peptide.
    """
    peptide = ""
    i = 0

    while i < len(sequence) and sequence[i] != '*':
        peptide += sequence[i]
        i += 1

    if peptide:
        return (peptide, 1, len(peptide))
    else:
        return None

def select_peptides(sequence):
    """
    Extracts peptides from a translated sequence, identifying sequences between stop codons (*).

    Args:
        sequence (str): Protein sequence.

    Returns:
        list of tuples: Each tuple contains (peptide, start, end) where:
            - peptide (str): The peptide sequence.
            - start (int): The starting position of the peptide in the sequence (1-indexed).
            - end (int): The ending position of the peptide.
    """
    peptides = []
    current_peptide = ""
    start = 0
    i = 0

    # Extract the peptide before the first stop codon
    while i < len(sequence) and sequence[i] != '*':
        current_peptide += sequence[i]
        i += 1

    if current_peptide:
        peptides.append((current_peptide, 1, len(current_peptide)))
        current_peptide = ""
    i += 1

    # Process peptides after the first stop codon
    while i < len(sequence):
        if sequence[i] == 'M':  # Start of a new peptide (initiator methionine)
            start = i
            current_peptide = 'M'
            i += 1
            while i < len(sequence) and sequence[i] != '*':
                current_peptide += sequence[i]
                i += 1
            if current_peptide:
                peptides.append((current_peptide, start + 1, i))
                current_peptide = ""
        i += 1

    return peptides

def process_translation(file_path, contig_sequences):
    """
    Processes a file containing translated peptide sequences and associates them with their respective frames.

    Args:
        file_path (str): Path to the file containing translated sequences.
        contig_sequences (dict): Dictionary of contig sequences.

    Returns:
        dict: A dictionary where keys are frame identifiers and values are lists of peptides with their metadata.
    """
    frames = {}
    current_frame = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_frame = line
                frames[current_frame] = []
            else:
                peptides = select_peptides(line)
                contig_id = current_frame.split("_frame")[0].replace(">", "")
                frames[current_frame].extend([(pep, start, end, contig_id) for (pep, start, end) in peptides])

    return frames

def map_to_contig_positions(peptides, contig_sequences):
    """
    Maps peptides to their original contig positions.

    Args:
        peptides (dict): Dictionary of peptides extracted from frames.
        contig_sequences (dict): Dictionary of contig sequences.

    Returns:
        dict: A dictionary of mapped peptides with their original contig positions and nucleotide sequences.
    """
    mapped_peptides = {}
    for frame, peptide_data in peptides.items():
        mapped_peptides[frame] = []
        for peptide, start, end, contig in peptide_data:
            contig_seq = contig_sequences.get(contig, "")
            if contig_seq:
                # Calculate nucleotide positions
                start_contig = 1 if start == 1 else start * 3
                end_contig = end * 3
                
                # Extract nucleotide sequence corresponding to the peptide
                sequence = contig_seq[start_contig-1:end_contig]

                # Ensure the mapped sequence is at least 25 nucleotides long
                if len(sequence) < 25:
                    remaining_length = 25 - len(sequence)
                    extra_nucleotides = contig_seq[end_contig:end_contig + remaining_length]
                    sequence += extra_nucleotides

                mapped_peptides[frame].append((peptide, start_contig, end_contig, sequence, contig))
    return mapped_peptides

def process_peptides_before_stop(file_path, contig_sequences):
    """
    Processes a file to extract only peptides before the first stop codon and map them to contigs.

    Args:
        file_path (str): Path to the file containing translated sequences.
        contig_sequences (dict): Dictionary of contig sequences.

    Returns:
        dict: A dictionary where keys are frame identifiers and values are peptides before the first stop codon.
    """
    frames = {}
    current_frame = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_frame = line
                frames[current_frame] = []
            else:
                # Extract the peptide before the first stop codon
                peptide_data = select_peptide_before_first_stop(line)
                if peptide_data:
                    peptide, start, end = peptide_data
                    contig_id = current_frame.split("_pep")[0].replace(">", "")
                    frames[current_frame].append((peptide, start, end, contig_id))

    return frames

def map_to_contig_positions_before_stop(peptides, contig_sequences):
    """
    Maps peptides before the first stop codon to their original contig positions.

    Args:
        peptides (dict): Dictionary of peptides extracted from frames.
        contig_sequences (dict): Dictionary of contig sequences.

    Returns:
        dict: A dictionary of mapped peptides with their original contig positions and nucleotide sequences.
    """
    mapped_peptides = {}
    for frame, peptide_data in peptides.items():
        mapped_peptides[frame] = []
        for peptide, start, end, contig in peptide_data:
            contig_seq = contig_sequences.get(contig, "")
            if contig_seq:
                # Calculate nucleotide positions
                start_contig = 1 if start == 1 else start * 3
                end_contig = end * 3
                
                # Extract nucleotide sequence corresponding to the peptide
                sequence = contig_seq[start_contig-1:end_contig]

                # Ensure the mapped sequence is at least 25 nucleotides long
                if len(sequence) < 25:
                    remaining_length = 25 - len(sequence)
                    extra_nucleotides = contig_seq[end_contig:end_contig + remaining_length]
                    sequence += extra_nucleotides

                mapped_peptides[frame].append((peptide, start_contig, end_contig, sequence, contig))
    return mapped_peptides

def load_contig_sequences(file_path):
    """
    Loads contig sequences from a FASTA file.

    Args:
        file_path (str): Path to the contig FASTA file.

    Returns:
        dict: A dictionary where keys are contig IDs and values are the corresponding sequences.
    """
    contig_sequences = {}
    current_contig = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_contig = line[1:]  # Remove the '>'
                contig_sequences[current_contig] = ""
            else:
                contig_sequences[current_contig] += line

    return contig_sequences

def main():
    """
    Main function to process the input files and map peptides to their respective contigs.
    The results are saved to an output file specified by the user.
    """
    if len(sys.argv) != 5:
        print("Usage: python script.py <translation_file> <contig_file> <output_file> <mode>")
        print("<mode>: 'before_stop' or 'all'")
        sys.exit(1)

    # Input file paths and output file
    translation_file = sys.argv[1]
    contig_file = sys.argv[2]
    output_file = sys.argv[3]
    mode = sys.argv[4]

    # Load contig sequences
    contig_sequences = load_contig_sequences(contig_file)

    if mode == 'before_stop':
        # Process peptides before the first stop codon
        selected_peptides = process_peptides_before_stop(translation_file, contig_sequences)
        print(selected_peptides)
        mapped_peptides = map_to_contig_positions(selected_peptides, contig_sequences)
        #print(mapped_peptides)
    elif mode == 'all':
        # Process all peptides
        selected_peptides = process_translation(translation_file, contig_sequences)
        print(selected_peptides)
        mapped_peptides = map_to_contig_positions(selected_peptides, contig_sequences)
        #print(mapped_peptides)
    else:
        print("Invalid mode. Use 'before_stop' or 'all'.")
        sys.exit(1)

    # Write results to the output file
    with open(output_file, "w") as output:
        output.write("Frame\tPeptide\tStart-End\tSequence\tContig\n")
        for frame, peptides in mapped_peptides.items():
            for peptide, start, end, sequence, contig in peptides:
                output.write(f"{frame}\t{peptide}\t{start}-{end}\t{sequence}\t{contig}\n")

if __name__ == "__main__":
    main()
