import sys
import numpy as np
from scipy.stats import binom_test

def read_table(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    table = [line.strip().split('\t') for line in lines]
    return table

def determine_major_phase(table, significance_threshold=0.05):
    new_table = [table[0] + ["p_value"]]  # Add the header and new "p_value" column to the new table
    for row in table[1:]:  # Start from the second line
        phases = [int(row[i]) for i in range(2, 5)]
        n = sum(phases)  # Total number of observations
        major_phase_index = np.argmax(phases)  # Index of the major phase
        major_phase_proportion = phases[major_phase_index] / n  # Proportion of the major phase
        p_value = binom_test(major_phase_proportion * n, n, p=1/3)  # Binomial test

        row.append(f"{p_value:.4f}")  # Add the p-value at the end of the row
        new_table.append(row)

    headers = new_table[0]
    data = new_table[1:]
    data_sorted = sorted(data, key=lambda x: float(x[-1]))  # Sort by p-value
    new_table = [headers] + data_sorted
    return new_table

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py path_to_table.tsv")
        sys.exit(1)

    file_path = sys.argv[1]
    table = read_table(file_path)
    new_table = determine_major_phase(table)

    # Display the new table with header and p-value column
    for row in new_table:
        print('\t'.join(row))
