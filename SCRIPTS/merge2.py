import pandas as pd
import re

# Step 1: Load initial files
fich1 = pd.read_csv("/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/KBFP/protocole_KBFP/pancreas_test/ORFpredtestpancreas/KmersFromContigsQuerySumPhaseSeqPvalueRSState", sep="\t", header=0, low_memory=False)
fich2 = pd.read_csv("/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/KBFP/protocole_KBFP/pancreas_test/all_contigsOfpeptides_RS+", sep="\t", header=None, names=["Frame", "Peptide", "Start-End", "Sequence", "Contig"])

# Step 2: Build a dictionary to manage correspondences
fich2_dict = {}
for _, row in fich2.iterrows():
    contig = row["Contig"]
    if contig not in fich2_dict:
        fich2_dict[contig] = []
    fich2_dict[contig].append(row["Sequence"])  # Add other columns if needed

# Step 3: Add a correspondence column to fich1
def get_correspondance(contig):
    if contig in fich2_dict:
        return ", ".join(fich2_dict[contig])  # Concatenate multiple correspondences
    return None

fich1["Correspondance"] = fich1["ID_contig"].map(get_correspondance)

# Step 4: Save the updated file
fich1.to_csv("fich1_updated.txt", sep="\t", index=False)
print("Update completed. File saved as fich1_updated.txt.")

# Step 5: Extract ID_contig and RSState columns
fich1_filtered = fich1[["ID_contig", "RSState"]]

# Step 6: Load files for contig correspondence
fich2_correspondance = pd.read_csv("/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/KBFP/protocole_KBFP/pancreas_test/all_contigsOfpeptides_RS+", sep="\t", header=None, names=["Frame", "Peptide", "Start-End", "Sequence", "Contig"])
fich2_correspondance_dict = fich2_correspondance.set_index("Peptide")["Contig"].to_dict()

# Step 7: Add the Contig column to fich1_filtered
fich1_filtered["Contig"] = fich1_filtered["ID_contig"].map(fich2_correspondance_dict)

# Step 8: Save the updated file
fich1_filtered.to_csv("fichier1_updated.txt", sep="\t", index=False)
print("Update completed. Result saved as fichier1_updated.txt.")

# Step 9: Load the updated file and the merging file
file1_path = "fichier1_updated.txt"
file2_path = "/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/DiscoveryPancreas_Patients_CellLines/final_annotated_pancreas_lensup35.tsv"

df1 = pd.read_csv(file1_path, sep="\t")
df2 = pd.read_csv(file2_path, sep="\t")

# Step 10: Extract elements before "_" from the 'Contig' column
df1['Extracted_Contig'] = df1['Contig'].apply(lambda x: re.split(r'_', str(x))[0])

# Step 11: Merge files
merged_df = df1.merge(df2, left_on='Extracted_Contig', right_on='tag', how='left')
merged_df.fillna('NA', inplace=True)

# Step 12: Rename columns and remove 'Extracted_Contig' column
merged_df.rename(columns={"ID_contig": "peptide", "Contig": "ID"}, inplace=True)
merged_df.drop(columns=["Extracted_Contig"], inplace=True)

# Step 13: Save the final merged file
merged_df.to_csv("merged_fichier_mod.txt", sep="\t", index=False)
print("Merging completed. Result saved as merged_fichier_2.txt.")

