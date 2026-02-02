import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_customized_curve(group, output_file_combined, marker_size=10, bar_width=0.3, group_spacing=1):
    plt.figure(figsize=(10, 5))

    # Color-blind-friendly colors for P1, P2, P3
    colors = ['#0072B2', '#E69F00', '#009E73']

    n = len(group)
    for i in range(n):
        tag = i // 3  # Each group of P1, P2, P3 gets the same tag
        position = tag * (3 * bar_width + group_spacing) + (i % 3) * bar_width
        sum_value = group['sum'].iloc[i]
        annotation = f"P{(i % 3) + 1}"

        # Use the corresponding color for each P1, P2, P3
        plt.bar(position, sum_value, width=bar_width, color=colors[i % 3], label=annotation if i < 3 else "")
        plt.text(position, sum_value, annotation, ha='center', va='bottom', fontsize=12)

    # Customize the plot
    plt.title(name, fontsize=22, fontweight='bold')
    plt.xlabel("Tag", fontweight='bold', fontsize=22)
    plt.ylabel("Count / Sum", fontweight='bold', fontsize=22)

    # Remove x-axis labels
    plt.xticks([])

    plt.yticks(fontsize=12, fontweight='bold')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)
    plt.tight_layout()

    # Save the combined plot as a PNG
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
            output_file_combined = f"{output_directory}/{name}_plot_hist.png"
            plot_customized_curve(group, output_file_combined, marker_size=4)
