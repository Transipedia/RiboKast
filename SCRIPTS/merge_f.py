import pandas as pd
import re
import argparse

# ==== STEP 1: PARSE ARGUMENTS ====
parser = argparse.ArgumentParser(description="Process and merge peptide and annotation data.")
parser.add_argument("-f1", "--fich1", required=True, help="Path to KmersFromContigsQuerySumPhaseSeqTranslatedPvalueRSState file")
parser.add_argument("-c", "--contigs", required=True, help="Path to contigsOfPeptides file (includes both RS+ and all_contigsOfpep)")
parser.add_argument("-a", "--annotation", required=True, help="Path to merged_annotation file")
parser.add_argument("-rsm", "--rsminus", required=True, help="Path to RS- file to include")
parser.add_argument("-o", "--output", required=True, help="Path to save the final merged output file")

args = parser.parse_args()

# ==== STEP 2: LOAD FILES ====
fich1 = pd.read_csv(args.fich1, sep="\t", header=0, low_memory=False)
fich2 = pd.read_csv(args.contigs, sep="\t", header=None, names=["Frame", "Peptide", "Start-End", "Sequence", "Contig"])
df2 = pd.read_csv(args.annotation, sep="\t")
rsminus_df = pd.read_csv(args.rsminus, sep="\t")

# ==== STEP 3: PREPARE RS+ DATA ====
fich2_dict = fich2.groupby("Contig")["Sequence"].apply(lambda x: ", ".join(x)).to_dict()
peptide_to_contig = fich2.set_index("Peptide")["Contig"].to_dict()

fich1["Correspondance"] = fich1["peptide"].map(fich2_dict)
fich1_filtered = fich1[["peptide", "RSState"]].copy()
fich1_filtered["Contig"] = fich1_filtered["peptide"].map(peptide_to_contig)
fich1_filtered['Extracted_Contig'] = fich1_filtered['Contig'].astype(str).apply(lambda x: re.split(r'_', x)[0])
print(fich1_filtered)
rsplus_merged = fich1_filtered.merge(df2, left_on='Extracted_Contig', right_on='tag', how='left').fillna('NA')
rsplus_merged.rename(columns={"ID_contig": "peptide", "Contig": "ID"}, inplace=True)
rsplus_merged.drop(columns=["Extracted_Contig"], inplace=True)

# ==== STEP 4: PREPARE RS- DATA ====
# Faire la jointure sur 'tag' sans transformer 'contig'
rsminus_merged = rsminus_df.merge(df2, on='tag', how='left').fillna('NA')

# Restaurer la colonne 'contig' de rsminus_df (elle peut être écrasée lors du merge)
if 'contig' in rsminus_df.columns:
    rsminus_merged['contig'] = rsminus_df['contig']

# Ajouter les colonnes manquantes si besoin
expected_columns = rsplus_merged.columns.tolist()
existing_columns = [col for col in expected_columns if col in rsminus_merged.columns]
missing_columns = [col for col in expected_columns if col not in rsminus_merged.columns]
for col in missing_columns:
    rsminus_merged[col] = 'NA'
rsminus_merged = rsminus_merged[expected_columns]

# ==== STEP 5: CONCAT RS+ AND RS- ====
final_df = pd.concat([rsplus_merged, rsminus_merged], ignore_index=True)

# ==== STEP 6: OUTPUT ====
final_df.to_csv(args.output, sep="\t", index=False)
print(f"Processing complete. Merged output saved to {args.output}")
