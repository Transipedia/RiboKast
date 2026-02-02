#!/bin/bash
#SBATCH --job-name="RiboKast"
#SBATCH --partition=ssfa -t 100:00:00 --mem 64G
#SBATCH --cpus-per-task=6

# ==== ACTIVATE CONDA ENVIRONMENT ====
source /home/safa.maddouri/miniconda3/bin/activate ribokast

# ==== LOAD CONFIG ====
source "/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_test/config.sh"

# ==== EXECUTE SCRIPTS ====

# Run with translation (-contig mode)
$BASE_DIR/run_RiboKast.sh -contig "$INDEX_DIR" "$SCRIPTS_DIR" "$CONTIGS_DIR" "$SIF_FILE" "$FASTA_FILE" "$PHASE_SHIFT" "$KMER_LEN"


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
