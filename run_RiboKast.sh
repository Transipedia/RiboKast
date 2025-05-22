#!/bin/bash
#SBATCH --job-name="kbfp"
#SBATCH --partition=ssfa -t 100:00:00 --mem 300G
#SBATCH --cpus-per-task=6

# ==== USAGE CHECK ====
if [ "$#" -lt 7 ]; then
    echo "Usage: $0 -contig|-orf <INDEX_DIR> <SCRIPTS_DIR> <FILES_DIR> <SIF_FILE> <FASTA_FILE> [<PHASE_SHIFT>] <KMER_LENGTH>"
    exit 1
fi

# ==== ARGUMENTS ====
MODE=$1            # -contig or -orf
INDEX_DIR=$2
SCRIPTS_DIR=$3
FILES_DIR=$4
SIF_FILE=$5
FASTA_FILE=$6

# Default phase shift to +1 if not provided
if [ "$#" -eq 7 ]; then
    PHASE_SHIFT="0"
    KMER_LENGTH=$7
else
    PHASE_SHIFT=$7
    KMER_LENGTH=$8
fi

# ==== VALIDATE TRANSLATION MODE ====
TRANSLATE=true  # Default: translation enabled

if [ "$MODE" = "-orf" ]; then
    TRANSLATE=false
elif [ "$MODE" != "-contig" ]; then
    echo "Invalid mode: use '-contig' to enable translation or '-orf' to disable it."
    exit 1
fi

# ==== INITIALIZATION ====
cd ${PBS_O_WORKDIR}

mkdir -p "$FILES_DIR"

# ==== KaMRat QUERY ====
apptainer exec --bind '/store:/store' -B '/data:/data' "$SIF_FILE" \
    kamrat query -idxdir "$INDEX_DIR" -fasta "$FASTA_FILE" -toquery median -withabsent -outpath "$FILES_DIR/out"

# ==== PROCESS OUTPUT FILE ====
python3 "$SCRIPTS_DIR/addid.py" "$FASTA_FILE" "$FILES_DIR/out" "$FILES_DIR/out_id"

# Retrieve header
HEADER=$(head -n 1 "$FILES_DIR/out_id")

# ==== FILTER RS+ / RS- ====
{
    echo -e "$HEADER" > "$FILES_DIR/RS+"
    awk -F'\t' 'NR > 1 {
        sum=0
        for(i=3; i<=NF; i++) sum+=$i
        if(sum != 0) print
    }' "$FILES_DIR/out_id" >> "$FILES_DIR/RS+"

    echo -e "$HEADER" > "$FILES_DIR/RS-"
    awk -F'\t' 'NR > 1 {
        sum=0
        for(i=3; i<=NF; i++) sum+=$i
        if(sum == 0) print
    }' "$FILES_DIR/out_id" >> "$FILES_DIR/RS-"
}

# Create RS+ / RS- FASTA files
awk 'NR > 1 {print ">"$1"\n"$2}' "$FILES_DIR/RS+" > "$FILES_DIR/RS+.fa"
awk 'NR > 1 {print ">"$1"\n"$2}' "$FILES_DIR/RS-" > "$FILES_DIR/RS-.fa"

# ==== GENERATE KMERS ====
python3 "$SCRIPTS_DIR/generate_kmers_fromFasta.py" "$FILES_DIR/RS+.fa" "$FILES_DIR/kmersFromContigs.fa" "$KMER_LENGTH"

# ==== SECOND KaMRat QUERY ====
apptainer exec --bind /store/EQUIPES/SSFA/MEMBERS/safa.maddouri "$SIF_FILE" \
    kamrat query -idxdir "$INDEX_DIR" -fasta "$FILES_DIR/kmersFromContigs.fa" -toquery median -withabsent -outpath "$FILES_DIR/KmersFromContigsQuery"

# ==== PHASING PREDICTION ====
python3 "$SCRIPTS_DIR/add_id_sum.py" "$FILES_DIR/kmersFromContigs.fa" "$FILES_DIR/KmersFromContigsQuery" "$FILES_DIR/KmersFromContigsQuerySum"
python3 "$SCRIPTS_DIR/phaseCount.py" "$FILES_DIR/KmersFromContigsQuerySum" "$FILES_DIR/KmersFromContigsQuerySumPhase" "$PHASE_SHIFT"
python3 "$SCRIPTS_DIR/add_colContFromFastaFile_arg.py" "$FILES_DIR/RS+.fa" "$FILES_DIR/KmersFromContigsQuerySumPhase" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeq"

# ==== TRANSLATION (Only if mode is '-contig') ====
if [ "$TRANSLATE" = true ]; then
    "$SCRIPTS_DIR/translate_st.sh" \
        "$FILES_DIR/KmersFromContigsQuerySumPhaseSeq" \
        "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqTranslated" \
        "$FILES_DIR/temp_fasta" \
        "$FILES_DIR/temp_result"
    rm "$FILES_DIR/temp_result" "$FILES_DIR/temp_fasta"

    # Binomial test after translation
    python3 "$SCRIPTS_DIR/binom_test.py" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqTranslated" > "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqTranslatedPvalue"
    python3 "$SCRIPTS_DIR/addRSState.py" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqTranslatedPvalue" "$FILES_DIR/RS-" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqTranslatedPvalueRSState"
else
    # Binomial test without translation
    python3 "$SCRIPTS_DIR/binom_test.py" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeq" > "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqPvalue"
    python3 "$SCRIPTS_DIR/addRSState.py" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqPvalue" "$FILES_DIR/RS-" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqPvalueRSState"

    # ==== HEADER UPDATE FOR -orf MODE ====
    TMP_HEADER_FILE="$FILES_DIR/tmp_header_replaced"
    awk 'NR==1 {
        printf "contigofpeptide\tpeptide";
        for(i=3; i<=NF; i++) {
            printf "\t%s", $i;
        }
        printf "\n";
        next;
    }
    { print }' "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqPvalueRSState" > "$TMP_HEADER_FILE"
    mv "$TMP_HEADER_FILE" "$FILES_DIR/KmersFromContigsQuerySumPhaseSeqPvalueRSState"
fi

echo "Pipeline execution completed successfully."
