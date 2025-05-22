#!/bin/bash

# Activate conda environment
source /home/safa.maddouri/miniconda3/bin/activate

# Load configuration
source "$(dirname "$0")/config.sh"

# Define working directory
FILES_DIR="$CONTIGS_DIR"
PLOTS_DIR="$FILES_DIR/plots"

# Create output directory if it doesn't exist
mkdir -p "$PLOTS_DIR"

# Execute post-processing scripts
python3 "$SCRIPTS_DIR/plotPhase.py" "$RSPLUS_KMERS" "$PLOTS_DIR"
python3 "$SCRIPTS_DIR/plotdist.py" "$FILES_DIR/KmersFromContigsQuerySum" "$PLOTS_DIR"
python3 "$SCRIPTS_DIR/merge_cont.py" "$PLOTS_DIR"
python3 "$SCRIPTS_DIR/plot_dis_phase_histogrames.py" "$FILES_DIR/KmersFromContigsQuerySum" "$PLOTS_DIR"

# Clean up generated intermediate plot files
find "$PLOTS_DIR" -type f \( -name "*_plot.png" -o -name "*_phase.png" \) -delete
