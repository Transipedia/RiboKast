import pandas as pd
import sys

def filter_contigs(input_file, output_file):
    # Charger le DataFrame à partir du fichier d'entrée
    df = pd.read_csv(input_file, sep='\t')
    
    # Groupement des lignes par contig
    #contig_groups = df.groupby(df['id'].str.extract(r'(pep_\d+)', expand=False))
    contig_groups = df.groupby(df['id'].str.extract(r'(.+)_kmer_\d+', expand=False))

    #contig_groups = df.groupby(df['id'].str.extract(r'^(\w+_\d+)', expand=False))

    # Liste pour stocker les contigs à garder
    contigs_to_keep = []
    
    # Parcourir les groupes de contigs
    for contig, group in contig_groups:
        # Calculer le nombre de lignes dans chaque groupe
        num_rows = len(group)
        # Calculer le nombre de lignes représentant au moins 60% du total
        #min_rows = 0.6 * num_rows
        #nbre de kmers dont le compte est différent de est supérieur ou égal à 21 pour chaque contig
        min_rows = 21
        # Vérifier si le nombre de lignes représentant au moins 60% est atteint
        print(group)
        if sum(group['sum'] > 0) >= min_rows:
            print(sum(group['sum'] > 0))
            contigs_to_keep.append(contig)
    
    # Filtrer le DataFrame original pour ne garder que les contigs sélectionnés
    #filtered_df = df[df['id'].str.extract(r'(pep_\d+)', expand=False).isin(contigs_to_keep)]
    filtered_df = df[df['id'].str.extract(r'(.+)_kmer_\d+', expand=False).isin(contigs_to_keep)]
    #filtered_df = df[df['id'].str.extract(r'^(\w+_\d+)', expand=False).isin(contigs_to_keep)]

    # Écrire le DataFrame filtré dans un fichier de sortie
    filtered_df.to_csv(output_file, sep='\t', index=False)

# Vérifier les arguments en ligne de commande
if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

# Récupérer les noms de fichiers à partir des arguments de ligne de commande
input_file = sys.argv[1]
output_file = sys.argv[2]

# Appeler la fonction pour filtrer les contigs
filter_contigs(input_file, output_file)

