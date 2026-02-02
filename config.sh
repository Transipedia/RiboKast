#!/bin/bash

# ==== BASE PATHS ====
BASE_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_test"
#OUT_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_test/RESULTS_CONTIGS_TEST_MELANOMA"
OUT_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_test_HKG/RESULTS_CDS_HKG"
INDEX_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_Indexes/Index_124CC_80%/RESULTS/kmerCount/Kamrat/index"
SCRIPTS_DIR="$BASE_DIR/SCRIPTS"
SIF_FILE="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/KaMRaT.sif"
FASTA_FILE="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_test_HKG/sequences_cds.fasta"
#FASTA_FILE="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKast_Indexes/Index_124CC_80%/RESULTS/kmerCount/Kamrat/merged-res-lensup35_test.fa"
#ANNOTATION_FILE="/datas/safa/melanoma/new_run_melanoma/annotation/merged_annotation.tsv"  # Leave commented if not used

# ==== KaMRaT PARAMETERS ====
KMER_LEN=25
KAMRAT_TOQUERY="mean"     # mean or median
KAMRAT_COUNTS="float"     # int or float
KAMRAT_WITHABSENT="1"     # 1 or 0

# ==== PHASE SHIFT ====
# 0 = no shift, 1 = +1, 2 = +2, etc.
PHASE_SHIFT="0"

# ==== PEPTIDE LENGTH FILTER (optional) ====
# If set, remove rows where peptide length is <= this value.
# Leave empty ("") to disable the filter.
PEPTIDE_LEN_MIN="9"

# ==== DIRECTORIES ==== DERIVED DIRECTORIES === (no need to change)

CONTIGS_DIR="$OUT_DIR/RESULTS_CONTIGS"
ORFPRED_DIR="$CONTIGS_DIR/RESULTS_ORFs"

# ==== FILES ====

FASTA_FILE_ORF="$CONTIGS_DIR/All_contig_of_peptides.fa"
OUTPUT_FILE="$ORFPRED_DIR/allORFsAnnotated.tsv"

# ==== MERGE INPUTS ====
RSPLUS_KMERS="$ORFPRED_DIR/KmersFromContigsQuerySumPhaseSeqPvalueRSState"
RSPLUS_CONTIGS="$CONTIGS_DIR/all_contigsOfpeptides_RS+"
RSMINUS_pep_CONTIGS="$CONTIGS_DIR/contigsOfPeptides_RS-"
RSMINUS_CONTIGS="$CONTIGS_DIR/RS-.tsv"
MERGE_OUTPUT_NO_ANNOT="$ORFPRED_DIR/ORFs_RSState.tsv"
