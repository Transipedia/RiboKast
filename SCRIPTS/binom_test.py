import sys
import numpy as np
from scipy.stats import binomtest

def read_table(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    table = [line.rstrip("\n").split('\t') for line in lines]
    return table

def determine_major_phase(table, significance_threshold=0.05):
    new_table = [table[0] + ["p_value"]]  # header + p_value column

    for row in table[1:]:
        phases = [int(row[i]) for i in range(2, 5)]
        n = sum(phases)

        # Avoid division by zero if a row has 0 total counts
        if n == 0:
            row.append("1.0000")
            new_table.append(row)
            continue

        major_phase_index = int(np.argmax(phases))
        k = int(phases[major_phase_index])  # successes = major phase count

        # Binomial test: H0 p=1/3
        p_value = binomtest(k=k, n=n, p=1/3).pvalue

        row.append(f"{p_value:.4f}")
        new_table.append(row)

    headers = new_table[0]
    data = new_table[1:]
    data_sorted = sorted(data, key=lambda x: float(x[-1]))  # sort by p_value
    return [headers] + data_sorted

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py path_to_table.tsv")
        sys.exit(1)

    file_path = sys.argv[1]
    table = read_table(file_path)
    new_table = determine_major_phase(table)

    for row in new_table:
        print('\t'.join(row))
