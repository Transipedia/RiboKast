import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_customized_curve(group, output_file_combined, initial_offset=40, offset_increment=10, sum_offset_increment=10, marker_size=10):
    plt.figure(figsize=(20, 8))  # Increase figure size

    # Plot count curves for each sample with vertical offset
    offset = initial_offset
    for i, column in enumerate(group.columns[2:-1]):  # Exclude 'tag' and 'Sum' columns
        plt.plot(group['tag'], group[column] + offset, marker='o', markersize=marker_size, linestyle='-', label=column, linewidth=2)
        offset += offset_increment

    # Plot sum curve with annotations P0, P1, P2
    sum_offset = initial_offset + (len(group.columns[1:-1]) - 1) * offset_increment + sum_offset_increment
    plt.plot(group['tag'], group['sum'] + sum_offset, marker='o', markersize=marker_size, linestyle='-', label='Sum', linewidth=2)
    for i in range(len(group['tag'])):
        tag = group['tag'].iloc[i]
        sum_value = group['sum'].iloc[i]
        annotation = f"P{(i % 3) + 1}"
        plt.text(tag, sum_value + sum_offset, annotation, ha='center', va='bottom', fontsize=16)

    # Customize the plot
    plt.title(name, fontsize=22, fontweight='bold')  # Increase title size
    plt.xlabel("Tag", fontweight='bold', fontsize=22)  # Increase x-axis label size
    plt.ylabel("Count / Sum", fontweight='bold', fontsize=22)  # Increase y-axis label size
    plt.xticks([])
    plt.yticks(fontsize=12, fontweight='bold')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
    plt.tight_layout()

    # Save the plot as a PNG file
    plt.savefig(output_file_combined)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_directory")
    else:
        input_file = sys.argv[1]
        output_directory = sys.argv[2]

        df = pd.read_csv(input_file, sep='\t')
        contig_groups = df.groupby(df['id'].str.extract(r'(.+)_kmer_\d+', expand=False))

        for name, group in contig_groups:
            output_file_combined = f"{output_directory}/{name}_plot.png"
            plot_customized_curve(group, output_file_combined, initial_offset=40, offset_increment=70, sum_offset_increment=20, marker_size=4)
