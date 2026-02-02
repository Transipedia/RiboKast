import pandas as pd
import sys

def process_contig_kmers(df, shift_value):
    results = {}  # Dictionary to store results by contig
    phases = ['p1', 'p2', 'p3']
    contig_groups = df.groupby(df['id'].str.extract(r'(.+)_kmer_\d+', expand=False))

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

        # Calculate the result based on the dominant phase and shift_value
        if shift_value == '+1':
            if dominant_phase == 'p1':
                result = 1 + 1
            elif dominant_phase == 'p2':
                result = 2 + 1
            elif dominant_phase == 'p3':
                result = 3 - 2
        elif shift_value == '0':
            if dominant_phase == 'p1':
                result = 1
            elif dominant_phase == 'p2':
                result = 2
            elif dominant_phase == 'p3':
                result = 3
        elif shift_value == '-1':
            if dominant_phase == 'p1':
                result = 1 + 2
            elif dominant_phase == 'p2':
                result = 2 - 1
            elif dominant_phase == 'p3':
                result = 3 - 1

        # Store the results for this contig
        results[contig] = {'P1': contig_kmer_counts['p1'], 'P2': contig_kmer_counts['p2'], 'P3': contig_kmer_counts['p3'],
                           'Dominant_Phase': dominant_phase, 'Functional_dominant_phase': result}  # Change here

    return results

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python script.py input_file output_file [shift_value]")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        shift_value = sys.argv[3] if len(sys.argv) == 4 else '0'  # Default shift_value = '0'

        df = pd.read_csv(input_file, sep='\t')

        # Process the kmers by contig
        results = process_contig_kmers(df, shift_value)
        print(results)
        # Write the results to the output file
        with open(output_file, 'w') as out_file:
            out_file.write("ID_contig\tP1\tP2\tP3\tDominant_Phase\tFunctional_dominant_phase\n")  # Change here
            for contig, result in results.items():
                out_file.write(f"{contig}\t{result['P1']}\t{result['P2']}\t{result['P3']}\t{result['Dominant_Phase']}\t{result['Functional_dominant_phase']}\n")
