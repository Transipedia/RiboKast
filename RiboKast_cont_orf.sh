#!/bin/bash
#SBATCH --job-name="RiboKast"
#SBATCH --partition=ssfa -t 100:00:00 --mem 64G
#SBATCH --cpus-per-task=6

# ==== ACTIVATE CONDA ENVIRONMENT ====
source /home/safa.maddouri/miniconda3/bin/activate ribokast

# ==== LOAD CONFIG ====
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
CONFIG_SH="${SUBMIT_DIR}/config.sh"
source "$CONFIG_SH"

# ==== EXECUTE SCRIPTS ====

# Run with translation (-contig mode)
$BASE_DIR/run_RiboKast.sh -contig "$INDEX_DIR" "$SCRIPTS_DIR" "$CONTIGS_DIR" "$SIF_FILE" "$FASTA_FILE" "$PHASE_SHIFT" "$KMER_LEN"

# Run ORF prediction
$BASE_DIR/ORFpred.sh "$CONTIGS_DIR" "$SCRIPTS_DIR"

# Run without translation (-orf mode)
$BASE_DIR/run_RiboKast.sh -orf "$INDEX_DIR" "$SCRIPTS_DIR" "$ORFPRED_DIR" "$SIF_FILE" "$FASTA_FILE_ORF" "$PHASE_SHIFT" "$KMER_LEN"

# ==== MERGE OUTPUT FILE ====
if [ -f "$ANNOTATION_FILE" ]; then
    echo "Annotation file found. Running annotated merge..."
    python3 "$SCRIPTS_DIR/merge_f.py" \
        -f1 "$RSPLUS_KMERS" \
        -c "$RSPLUS_CONTIGS" \
        -a "$ANNOTATION_FILE" \
        -rsm "$RSMINUS_CONTIGS" \
        -o "$OUTPUT_FILE"
else
    echo "Annotation file not found. Running merge without annotation..."
    python3 "$SCRIPTS_DIR/merge_no_ann.py" \
        "$RSPLUS_KMERS" \
        "$RSPLUS_CONTIGS" \
        "$RSMINUS_pep_CONTIGS" \
        "$MERGE_OUTPUT_NO_ANNOT"
fi



# ==== OPTIONAL: REMOVE LINES WHERE peptide LENGTH <= PEPTIDE_LEN_MIN ====

# ==== SELECT FINAL OUTPUT FILE (annotated or not) ====
FINAL_OUT=""
if [ -f "$OUTPUT_FILE" ]; then
    FINAL_OUT="$OUTPUT_FILE"
elif [ -f "$MERGE_OUTPUT_NO_ANNOT" ]; then
    FINAL_OUT="$MERGE_OUTPUT_NO_ANNOT"
fi

# ==== OPTIONAL: REMOVE LINES WHERE peptide LENGTH <= PEPTIDE_LEN_MIN ====
if [ -n "${PEPTIDE_LEN_MIN:-}" ]; then
    if [ -n "$FINAL_OUT" ] && [ -f "$FINAL_OUT" ]; then
        echo "Filtering: removing rows where peptide length <= $PEPTIDE_LEN_MIN (file: $FINAL_OUT)"

        tmp_out="${FINAL_OUT}.tmp_lenfilter"

        awk -F'\t' -v OFS='\t' -v minlen="$PEPTIDE_LEN_MIN" '
            NR==1{
                col=0
                for(i=1;i<=NF;i++){
                    h=tolower($i)
                    gsub(/\r/,"",h)
                    if(h=="peptide"){ col=i; break }
                }
                if(col==0){
                    print "WARNING: column \"peptide\" not found in header. No filtering applied." > "/dev/stderr"
                    print $0
                    next
                }
                peptide_col=col
                print $0
                next
            }
            {
                p=$peptide_col
                gsub(/\r/,"",p)
                # KEEP only peptides with length > minlen
                if(length(p) >= minlen) print $0
            }
        ' "$FINAL_OUT" > "$tmp_out"

        if [ -s "$tmp_out" ]; then
            mv "$tmp_out" "$FINAL_OUT"
            echo "Filtering done. Output updated: $FINAL_OUT"
        else
            echo "WARNING: filtering produced an empty file. Keeping original."
            rm -f "$tmp_out"
        fi
    else
        echo "WARNING: no final output file found (OUTPUT_FILE or MERGE_OUTPUT_NO_ANNOT). Skipping peptide-length filtering."
    fi
else
    echo "No PEPTIDE_LEN_MIN set => no peptide-length filtering."
fi


# ============================
# CLEAN TEMP / INTERMEDIATE FILES
# ============================

# ---- 1) Intermediates from run_RiboKast.sh (remove in BOTH contigs + orfpred dirs)
for d in "$CONTIGS_DIR" "$ORFPRED_DIR"; do
  rm -f \
    "$d/out" \
    "$d/RS+" \
    "$d/RS-" \
    "$d/kmersFromContigs.fa" \
    "$d/KmersFromContigsQuery" \
    "$d/KmersFromContigsQuerySumPhase" \
    "$d/KmersFromContigsQuerySumPhaseSeq" \
    "$d/KmersFromContigsQuerySumPhaseSeqTranslatedPvalue" \
    "$d/KmersFromContigsQuerySumPhaseSeqTranslated"
done

# ---- 2) ORFpred intermediates (these live in ORFPRED_DIR)
rm -f \
  "$CONTIGS_DIR/all_contigsOfpeptides_RS+" \
  "$CONTIGS_DIR/contigsOfPeptides_RS-" \
  "$CONTIGS_DIR/contigsOfPeptides_RS+P-" \
  "$CONTIGS_DIR/contigsOfPeptides_RS+P+" \
  "$CONTIGS_DIR/RS-_pep.fa" \
  "$CONTIGS_DIR/RS+P-_pep.fa" \
  "$CONTIGS_DIR/RS.fa-" \
  "$CONTIGS_DIR/RS+P+_pep.fa" \
  "$CONTIGS_DIR/RS+P-.fa"
