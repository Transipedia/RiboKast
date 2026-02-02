#!/bin/bash

# ==== CHECK USAGE ====
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <FILES_DIR> <SCRIPTS_DIR> <INPUT_FILE>"
    exit 1
fi

# ==== ARGUMENTS ====
FILES_DIR=$1  # Main directory path
SCRIPTS_DIR=$2  # Scripts directory
INPUT_FILE=$3  # Input file

# Output files (inside FILES_DIR)
RS_P_PLUS_FA="$FILES_DIR/RS+P+.fa"
RS_P_PLUS_PEP_FA="$FILES_DIR/RS+P+_pep.fa"
RS_P_MINUS_FA="$FILES_DIR/RS+P-.fa"
RS_P_MINUS_PEP_FA="$FILES_DIR/RS+P-_pep.fa"
FINAL_CONTIGS="$FILES_DIR/All_contig_of_peptides.fa"

# ==== GENERATE FASTA FILES USING AWK ====
awk -v rs_p_plus_fa="$RS_P_PLUS_FA" -v rs_p_plus_pep_fa="$RS_P_PLUS_PEP_FA" -v rs_p_minus_fa="$RS_P_MINUS_FA" '{
    if ($NF == "RS+P+") {
        # Create FASTA file for RS+P+ DNA sequences
        print ">" $2 "\n" $1 > rs_p_plus_fa
        # Create FASTA file for RS+P+ protein sequences
        print ">" $2 "_pep" "\n" $8 > rs_p_plus_pep_fa
    } else if ($NF == "RS+P-") {
        # Create FASTA file for RS+P- DNA sequences
        print ">" $2 "\n" $1 > rs_p_minus_fa
    }
}' "$INPUT_FILE"

# ==== TRANSLATE RS+P- DNA SEQUENCES TO PROTEIN ====
seqkit translate -F -f 1,2,3 "$RS_P_MINUS_FA" --line-width 7000 > "$RS_P_MINUS_PEP_FA"

# ==== EXTRACT ORFs ====
python3 "$SCRIPTS_DIR/getORF.py" "$RS_P_PLUS_PEP_FA" "$RS_P_PLUS_FA" "$FILES_DIR/contigsOfPeptides_RS+P+" before_stop
python3 "$SCRIPTS_DIR/getORF.py" "$RS_P_MINUS_PEP_FA" "$RS_P_MINUS_FA" "$FILES_DIR/contigsOfPeptides_RS+P-" all

# ==== FILTER CONTIGS BASED ON LENGTH ====
awk 'NR > 1 && length($2) >= 9 {print ">"$2"\n"$4}' "$FILES_DIR/contigsOfPeptides_RS+P+" > "$FILES_DIR/contigsOfpeptides_RS+P+.fa"
awk 'NR > 1 && length($2) >= 9 {print ">"$2"\n"$4}' "$FILES_DIR/contigsOfPeptides_RS+P-" > "$FILES_DIR/contigsOfpeptides_RS+P-.fa"

# ==== MERGE RS+P+ AND RS+P- CONTIG FILES INTO ONE ====
cat "$FILES_DIR/contigsOfpeptides_RS+P+.fa" "$FILES_DIR/contigsOfpeptides_RS+P-.fa" > "$FINAL_CONTIGS"

echo "Pipeline completed. Final file: $FINAL_CONTIGS"
