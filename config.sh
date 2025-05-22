#!/bin/bash

# ==== BASE PATHS ====
BASE_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/KBFP/protocole_KBFP/protocole_KBFP_20_02_2025/RiboKast"
INDEX_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/KBFP/index_RS/IndexKamrat_20mer"
SCRIPTS_DIR="$BASE_DIR/SCRIPTS"
SIF_FILE="/store/EQUIPES/SSFA/daniel.shared/Bin/KaMRaT.sif"

# ==== DIRECTORIES ====
CONTIGS_DIR="$BASE_DIR/pancreas_test"
ORFPRED_DIR="$CONTIGS_DIR/ORFpredtestpancreas"

# ==== FILES ====
FASTA_FILE="$BASE_DIR/test"
FASTA_FILE_ORF="$CONTIGS_DIR/All_contig_of_peptides.fa"
#ANNOTATION_FILE="$BASE_DIR/final_annotated_pancreas_lensup35.tsv"  # Leave commented if not used
OUTPUT_FILE="$ORFPRED_DIR/allORFsAnnotated.tsv"

# ==== MERGE INPUTS ====
RSPLUS_KMERS="$CONTIGS_DIR/KmersFromContigsQuerySumPhaseSeqTranslatedPvalueRSState"
RSPLUS_CONTIGS="$CONTIGS_DIR/all_contigsOfpeptides_RS+"
RSMINUS_CONTIGS="$CONTIGS_DIR/contigsOfPeptides_RS-"
MERGE_OUTPUT_NO_ANNOT="$ORFPRED_DIR/ORFs_RSState.tsv"
