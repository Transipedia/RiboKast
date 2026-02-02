#!/usr/bin/env bash
#SBATCH --job-name="RiboKast_postprocess"
#SBATCH --partition=ssfa
#SBATCH -t 100:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=6

set -euo pipefail
IFS=$'\n\t'

# ---- Conda + config ----
source /home/safa.maddouri/miniconda3/bin/activate ribokast
source "/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_test/config.sh"

RESULTS_DIR="$CONTIGS_DIR"
PLOTS_DIR="$RESULTS_DIR/plots"
mkdir -p "$PLOTS_DIR"

# ---- Inputs ----
INPUT_RSSTATE="$RESULTS_DIR/KmersFromContigsQuerySumPhaseSeqTranslatedPvalueRSState"
IN_TSV="$RESULTS_DIR/KmersFromContigsQuerySum"

# ---- Output FASTA we will generate here ----
FASTA_RS="$RESULTS_DIR/RS+P+.fa"

[[ -s "$INPUT_RSSTATE" ]] || { echo "[ERROR] Missing/empty input: $INPUT_RSSTATE" >&2; exit 1; }
[[ -s "$IN_TSV"        ]] || { echo "[ERROR] Missing/empty input: $IN_TSV" >&2; exit 1; }

# =========================
# 0) Generate RS+P+.fa from the RSState table
# =========================
# Expected columns (as in your README example):
#   $1 = contig sequence
#   $2 = ID_contig
#   ...
#   $NF = RSState (RS+P+, RS+P-, RS-)
#
# If your file has a header line, we skip it (NR==1) if it contains "RSState" or "ID_contig".
awk -v out_fa="$FASTA_RS" '
BEGIN {
  # overwrite output to avoid appending from previous runs
  print "" > out_fa
}
NR==1 {
  # Heuristic: if first line looks like a header, skip it
  if ($0 ~ /RSState/ || $0 ~ /ID_contig/ || $0 ~ /^contig[ \t]/) next
}
{
  if ($NF == "RS+P+") {
    # FASTA format: >ID then sequence
    print ">" $2 "\n" $1 >> out_fa
  }
}
END {
  # remove the initial empty line we wrote in BEGIN if file got content
}' "$INPUT_RSSTATE"

# Remove possible leading empty line (safe)
sed -i '1{/^$/d;}' "$FASTA_RS" || true

[[ -s "$FASTA_RS" ]] || { echo "[ERROR] RS+P+ FASTA is empty: $FASTA_RS" >&2; exit 1; }
echo "[INFO] Generated: $FASTA_RS (RS+P+ contigs)"

# =========================
# Temporary intermediate files (auto-removed)
# =========================
IDS_FILE="$(mktemp "$RESULTS_DIR/.RSplusPplus_ids.XXXXXX")"
FILTERED_TSV="$(mktemp "$RESULTS_DIR/.KmersFromContigsQuerySum_RSplusPplus.XXXXXX")"

cleanup() { rm -f "$IDS_FILE" "$FILTERED_TSV"; }
trap cleanup EXIT

# 1) Extract contig IDs from FASTA
grep '^>' "$FASTA_RS" | sed 's/^>//' | cut -d' ' -f1 > "$IDS_FILE"
[[ -s "$IDS_FILE" ]] || { echo "[ERROR] No IDs extracted from FASTA headers in: $FASTA_RS" >&2; exit 1; }

# 2) Filter while ALWAYS keeping the header from input
awk '
NR==FNR { keep[$1]=1; next }
FNR==1  { print $0; next }
{
  base=$1
  sub(/_kmer_[0-9]+$/, "", base)
  if (keep[base]) print $0
}
' "$IDS_FILE" "$IN_TSV" > "$FILTERED_TSV"

echo "[INFO] Using filtered temporary file: $FILTERED_TSV"
echo "[INFO] Header:"
head -n 1 "$FILTERED_TSV"

# 3) Plots on filtered file
python3 "$SCRIPTS_DIR/plotdist.py" "$FILTERED_TSV" "$PLOTS_DIR"
python3 "$SCRIPTS_DIR/plot_dis_phase_histogrames.py" "$FILTERED_TSV" "$PLOTS_DIR"
Rscript "$SCRIPTS_DIR/heatmaps.R" "$FILTERED_TSV" "$PLOTS_DIR"
echo "[INFO] Done. Intermediate files cleaned."

