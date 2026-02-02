#!/usr/bin/env python3
import sys

# ----------------------------
# Peptide selection functions
# ----------------------------

def select_peptide_before_first_stop(sequence: str):
    """
    Extract peptide sequence before the first stop codon (*).

    Returns:
        (peptide, start, end) with start/end 1-indexed AA positions, or None.
    """
    peptide = ""
    i = 0
    while i < len(sequence) and sequence[i] != '*':
        peptide += sequence[i]
        i += 1

    if peptide:
        return (peptide, 1, len(peptide))
    return None


def select_peptides(sequence: str):
    """
    Extract peptides from a translated sequence:
    - peptide before the first stop codon
    - then peptides starting at each 'M' after a stop, until next stop

    Returns:
        list of (peptide, start, end) with start/end 1-indexed AA positions.
    """
    peptides = []
    current_peptide = ""
    i = 0

    # peptide before first stop
    while i < len(sequence) and sequence[i] != '*':
        current_peptide += sequence[i]
        i += 1
    if current_peptide:
        peptides.append((current_peptide, 1, len(current_peptide)))
        current_peptide = ""

    # skip the first stop (if any)
    if i < len(sequence) and sequence[i] == '*':
        i += 1

    # peptides after first stop: start at each M
    while i < len(sequence):
        if sequence[i] == 'M':
            start_i = i  # 0-based index in AA string
            current_peptide = 'M'
            i += 1
            while i < len(sequence) and sequence[i] != '*':
                current_peptide += sequence[i]
                i += 1
            # end position is i (0-based stop index), so AA end = i (1-indexed)
            peptides.append((current_peptide, start_i + 1, i))
            current_peptide = ""
        i += 1

    return peptides


# ----------------------------
# FASTA loading
# ----------------------------

def load_contig_sequences(file_path: str):
    """
    Loads contig nucleotide sequences from FASTA.

    Keys are contig IDs (header without '>').
    """
    contig_sequences = {}
    current_contig = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_contig = line[1:].strip()
                contig_sequences[current_contig] = ""
            else:
                if current_contig is None:
                    continue
                contig_sequences[current_contig] += line

    return contig_sequences


# ----------------------------
# Helpers
# ----------------------------

def normalize_contig_id(header_line: str):
    """
    Normalize contig ID between translation headers and contig fasta headers.

    Works with headers like:
      >AAAA..._posSeq_pep
      >AAAA..._posSeq_frame1_pep
      >AAAA..._negSeq_frame2
    Produces base contig id by stripping:
      - leading '>'
      - any '_frame...' suffix
      - any '_pep' suffix
    """
    cid = header_line[1:] if header_line.startswith(">") else header_line
    cid = cid.split("_frame")[0]
    cid = cid.split("_pep")[0]
    return cid


# ----------------------------
# Translation processing
# ----------------------------

def process_translation(file_path: str):
    """
    Parse translated sequences FASTA-like file:
      >header
      AASEQ...

    Returns:
      dict frame_header -> list of (peptide, aa_start, aa_end, contig_id)
    """
    frames = {}
    current_frame = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_frame = line
                frames[current_frame] = []
            else:
                peptides = select_peptides(line)
                contig_id = normalize_contig_id(current_frame)
                frames[current_frame].extend(
                    [(pep, start, end, contig_id) for (pep, start, end) in peptides]
                )

    return frames


def process_peptides_before_stop(file_path: str):
    """
    Parse translated sequences and keep ONLY peptide before first stop.

    Returns:
      dict frame_header -> list of (peptide, aa_start, aa_end, contig_id)
    """
    frames = {}
    current_frame = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_frame = line
                frames[current_frame] = []
            else:
                peptide_data = select_peptide_before_first_stop(line)
                if peptide_data:
                    peptide, start, end = peptide_data
                    contig_id = normalize_contig_id(current_frame)
                    frames[current_frame].append((peptide, start, end, contig_id))

    return frames


# ----------------------------
# Mapping peptides to contig positions
# ----------------------------

def map_to_contig_positions(peptides_by_frame: dict, contig_sequences: dict, min_nt_len: int = 25, debug: bool = False):
    """
    Map AA positions to nucleotide positions in the original contig.
    Assumes frame already accounted for in translation, so AA->NT is *3 mapping.

    Returns:
      dict frame_header -> list of (peptide, nt_start, nt_end, nt_sequence, contig_id)
    """
    mapped = {}

    for frame, peptide_data in peptides_by_frame.items():
        mapped[frame] = []

        for peptide, aa_start, aa_end, contig in peptide_data:
            contig_seq = contig_sequences.get(contig, "")

            # Try fallback: sometimes contig fasta headers still contain suffixes
            # (rare) -> attempt to match by partial key
            if not contig_seq and debug:
                print(f"[DEBUG] contig not found: {contig} (frame {frame})")

            if not contig_seq:
                continue

            # Convert AA positions (1-indexed) to NT positions (1-indexed)
            nt_start = 1 if aa_start == 1 else aa_start * 3
            nt_end = aa_end * 3

            # Extract nucleotide sequence corresponding to peptide
            seq = contig_seq[nt_start - 1: nt_end]

            # Pad to minimum nucleotide length (e.g., 25) by extending downstream
            if len(seq) < min_nt_len:
                remaining = min_nt_len - len(seq)
                seq += contig_seq[nt_end: nt_end + remaining]

            mapped[frame].append((peptide, nt_start, nt_end, seq, contig))

    return mapped


# ----------------------------
# Main
# ----------------------------

def main():
    if len(sys.argv) != 5:
        print("Usage: python getORF.py <translation_file> <contig_fasta> <output_file> <mode>")
        print("<mode>: 'before_stop' or 'all'")
        sys.exit(1)

    translation_file = sys.argv[1]
    contig_file = sys.argv[2]
    output_file = sys.argv[3]
    mode = sys.argv[4]

    contig_sequences = load_contig_sequences(contig_file)

    if mode == "before_stop":
        selected = process_peptides_before_stop(translation_file)
        mapped = map_to_contig_positions(selected, contig_sequences)
    elif mode == "all":
        selected = process_translation(translation_file)
        mapped = map_to_contig_positions(selected, contig_sequences)
    else:
        print("Invalid mode. Use 'before_stop' or 'all'.")
        sys.exit(1)

    with open(output_file, "w") as out:
        out.write("Frame\tPeptide\tStart-End\tSequence\tContig\n")
        for frame, pep_list in mapped.items():
            for peptide, nt_start, nt_end, seq, contig in pep_list:
                out.write(f"{frame}\t{peptide}\t{nt_start}-{nt_end}\t{seq}\t{contig}\n")


if __name__ == "__main__":
    main()
