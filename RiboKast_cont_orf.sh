#!/bin/bash
#SBATCH --job-name="RiboKast"
#SBATCH --partition=ssfa -t 100:00:00 --mem 200G
#SBATCH --cpus-per-task=6

# ==== ACTIVATE CONDA ENVIRONMENT ====
source /home/safa.maddouri/miniconda3/bin/activate

# ==== LOAD CONFIG ====
source "$(dirname "$0")/config.sh"

# ==== EXECUTE SCRIPTS ====

# Run with translation (-contig mode)
$BASE_DIR/run_RiboKast.sh -contig "$INDEX_DIR" "$SCRIPTS_DIR" "$CONTIGS_DIR" "$SIF_FILE" "$FASTA_FILE" 20

# Run ORF prediction
$BASE_DIR/ORFpred.sh "$CONTIGS_DIR" "$SCRIPTS_DIR"

# Run without translation (-orf mode)
$BASE_DIR/run_RiboKast.sh -orf "$INDEX_DIR" "$SCRIPTS_DIR" "$ORFPRED_DIR" "$SIF_FILE" "$FASTA_FILE_ORF" 20

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
        "$RSMINUS_CONTIGS" \
        "$MERGE_OUTPUT_NO_ANNOT"
fi

echo "Pipeline execution completed successfully."
