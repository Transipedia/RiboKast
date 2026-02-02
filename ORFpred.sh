#!/bin/bash

# ==== CHECK USAGE ====
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <FILES_DIR> <SCRIPTS_DIR>"
    exit 1
fi

# ==== ARGUMENTS ====
FILES_DIR=$1  # Main directory path
SCRIPTS_DIR=$2  # Scripts directory

# ==== INPUT FILE ====
INPUT_FILE="$FILES_DIR/KmersFromContigsQuerySumPhaseSeqTranslatedPvalueRSState"

# ==== OUTPUT FILES ====
RS_P_PLUS_FA="$FILES_DIR/RS+P+.fa"
RS_P_PLUS_PEP_FA="$FILES_DIR/RS+P+_pep.fa"
RS_P_MINUS_FA="$FILES_DIR/RS+P-.fa"
RS_P_MINUS_PEP_FA="$FILES_DIR/RS+P-_pep.fa"
FINAL_CONTIGS="$FILES_DIR/All_contig_of_peptides.fa"
RS_MINUS="$FILES_DIR/RS-"
RS_MINUS_FA="$FILES_DIR/RS.fa-"
RS_MINUS_PEP_FA="$FILES_DIR/RS-_pep.fa"
RS_MINUS_TABLE="$FILES_DIR/RS-.tsv"

# ==== INITIALIZE OUTPUT FILES ====
: > "$RS_P_PLUS_FA"
: > "$RS_P_MINUS_FA"
: > "$RS_MINUS_FA"
: > "$RS_P_PLUS_PEP_FA"
: > "$RS_P_MINUS_PEP_FA"
: > "$RS_MINUS_PEP_FA"
# ==== GENERATE FASTA FILES USING AWK ====
awk -v rs_p_plus_fa="$RS_P_PLUS_FA" -v rs_p_plus_pep_fa="$RS_P_PLUS_PEP_FA" -v rs_p_minus_fa="$RS_P_MINUS_FA" '{
    if ($NF == "RS+P+") {
        print ">" $2 "\n" $1 > rs_p_plus_fa
        print ">" $2 "_pep" "\n" $8 > rs_p_plus_pep_fa
    } else if ($NF == "RS+P-") {
        print ">" $2 "\n" $1 > rs_p_minus_fa
    }
}' "$INPUT_FILE"

awk 'NR > 1 {print ">"$1"\n"$2}' "$RS_MINUS" > "$RS_MINUS_FA"

# ==== TRANSLATE RS+P-/RS- DNA SEQUENCES TO PROTEIN ====
seqkit translate -F -f 1,2,3 "$RS_P_MINUS_FA" --line-width 7000 > "$RS_P_MINUS_PEP_FA"
seqkit translate -F -f 1,2,3 "$RS_MINUS_FA" --line-width 7000 > "$RS_MINUS_PEP_FA"

# ==== EXTRACT ORFs ====
python3 "$SCRIPTS_DIR/getORF.py" "$RS_P_PLUS_PEP_FA" "$RS_P_PLUS_FA" "$FILES_DIR/contigsOfPeptides_RS+P+" all
python3 "$SCRIPTS_DIR/getORF.py" "$RS_P_MINUS_PEP_FA" "$RS_P_MINUS_FA" "$FILES_DIR/contigsOfPeptides_RS+P-" all
python3 "$SCRIPTS_DIR/getORF.py" "$RS_MINUS_PEP_FA" "$RS_MINUS_FA" "$FILES_DIR/contigsOfPeptides_RS-" all

cat "$FILES_DIR/contigsOfPeptides_RS+P+" "$FILES_DIR/contigsOfPeptides_RS+P-" > "$FILES_DIR/all_contigsOfpeptides_RS+"

# ==== TREAT RS- ====
awk '
BEGIN {FS=OFS="\t"}

FNR==1 && NR!=FNR {
    print "peptide", "ID", "contig", "RSState", "tag"
}
FNR==NR {
    if ($0 ~ /^>/) {
        id = substr($0, 2)
        next
    }
    contig_seq[id] = $0
    next
}
NR > FNR {
    if (FNR == 1) next
    peptide = $2
    contig_id = $5
    seq = contig_seq[contig_id]
    split(contig_id, parts, "_")
    tag = parts[1]
    print peptide, contig_id, seq, "RS-", tag
}
' "$RS_MINUS_FA" "$FILES_DIR/contigsOfPeptides_RS-" > "$RS_MINUS_TABLE"

# ==== FILTER CONTIGS BASED ON LENGTH ====
awk 'NR > 1 && length($2) >= 9 {print ">"$2"\n"$4}' "$FILES_DIR/contigsOfPeptides_RS+P+" > "$FILES_DIR/contigsOfpeptides_RS+P+.fa"
awk 'NR > 1 && length($2) >= 9 {print ">"$2"\n"$4}' "$FILES_DIR/contigsOfPeptides_RS+P-" > "$FILES_DIR/contigsOfpeptides_RS+P-.fa"
awk 'NR > 1 && length($2) >= 9 {print ">"$2"\n"$3}' "$RS_MINUS_TABLE" > "$FILES_DIR/contigsOfpeptides_RS-.fa"

# ==== MERGE RS+P+ AND RS+P- CONTIG FILES INTO ONE ====
cat "$FILES_DIR/contigsOfpeptides_RS+P+.fa" "$FILES_DIR/contigsOfpeptides_RS+P-.fa" > "$FINAL_CONTIGS"

echo "Pipeline completed. Final file: $FINAL_CONTIGS"
