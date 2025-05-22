import pandas as pd
import sys

def process_contig_kmers(df):
    results = {}  # Dictionary to store results by contig
    phases = ['p1', 'p2', 'p3']
    # Group the data by contig
    contig_groups = df.groupby(df['id'].str.extract(r'(pep_\d+)', expand=False))
    contig_groups = df.groupby(df['id'].str.extract(r'(.+)_kmer_\d+', expand=False))
    # Iterate over each group of kmers by contig
    for contig, group_df in contig_groups:
        # Initialize counters for this contig
        contig_kmer_counts = {phase: 0 for phase in phases}
        max_sum_counts = {phase: 0 for phase in phases}

        # Iterate over each group of 3 kmers
        for i in range(0, len(group_df), 3):
            # Extract the Sum values for each kmer in the group of 3 kmers
            sum_values = group_df.iloc[i:i+3]['sum'].values

            # Find the index of the maximum value
            max_index = sum_values.argmax() if any(sum_values) else None

            # Update the counters for the phase with the maximum value
            if max_index is not None:
                contig_kmer_counts[phases[max_index]] += 1
                max_sum_counts[phases[max_index]] += 1

        # Check if all counters are zero
        if all(value == 0 for value in contig_kmer_counts.values()):
            continue  # Move to the next contig if all counters are zero

        # Find the dominant phase
        dominant_phase = max(max_sum_counts, key=max_sum_counts.get)

        # Calculate the result based on the dominant phase
        if dominant_phase == 'p1':
            result = 2
        elif dominant_phase == 'p2':
            result = 3
        elif dominant_phase == 'p3':
            result = 1

        # Store the results for this contig
        results[contig] = {'P1': contig_kmer_counts['p1'], 'P2': contig_kmer_counts['p2'], 'P3': contig_kmer_counts['p3'],
                           'Dominant_Phase': dominant_phase, 'Phase+1': result}

    return results

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        df = pd.read_csv(input_file, sep='\t')
        # Process kmers by contig
        results = process_contig_kmers(df)
        print(results)
        # Write the results to the output file
        with open(output_file, 'w') as out_file:
            out_file.write("ID_contig\tP1\tP2\tP3\tDominant_Phase\tPhase+1\n")
            for contig, result in results.items():   
                out_file.write(f"{contig}\t{result['P1']}\t{result['P2']}\t{result['P3']}\t{result['Dominant_Phase']}\t{result['Phase+1']}\n")
