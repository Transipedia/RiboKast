import sys
import csv

# === Input files ===
rsplus_input = sys.argv[1]          # KmersFromContigsQuerySumPhaseSeqPvalueRSState
rsplus_mapping = sys.argv[2]        # all_contigsOfpeptides_RS+
rsminus_file = sys.argv[3]          # contigsOfPeptides_RS-
output_path = sys.argv[4]           # final output

# === Step 1: build peptide → contig_id mapping from RS+ ===
peptide_to_contig_rsplus = {}

with open(rsplus_mapping, newline='') as map_file:
    reader = csv.DictReader(map_file, delimiter='\t')
    for row in reader:
        peptide = row['Peptide'].strip()
        contig = row['Contig'].strip().lstrip('>')
        if peptide and contig:
            peptide_to_contig_rsplus[peptide] = contig

# === Step 2: process RS+ file and replace column 0 ===
data_rows = []

with open(rsplus_input, newline='') as infile:
    reader = csv.reader(infile, delimiter='\t')
    header = next(reader)
    col_count = len(header)
    header[0] = 'Contig_id'
    data_rows.append(header)

    for row in reader:
        peptide = row[1].strip()
        contig = peptide_to_contig_rsplus.get(peptide, 'UNKNOWN')
        row[0] = contig
        data_rows.append(row)

# === Step 3: process RS- ===
with open(rsminus_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        peptide = row['Peptide'].strip()
        contig = row['Contig'].strip().lstrip('>')
        new_row = ['NA'] * col_count
        new_row[0] = contig
        new_row[1] = peptide
        new_row[-1] = 'RS-'  # RSState = last column
        data_rows.append(new_row)

# === Step 4: write the final combined file ===
with open(output_path, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerows(data_rows)
