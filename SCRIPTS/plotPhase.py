import os
import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_histogram(row, output_dir):
    # Create the plot with a 3/2 aspect ratio
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.bar(['P1', 'P2', 'P3'], row[['P1', 'P2', 'P3']], color='skyblue', alpha=0.5)

    # Customize the plot
    ax.set_xlabel('Sum by phase', fontsize=22)
    ax.set_title(f'{row["ID_contig"]} - pvalue={row["p_value"]}', weight='bold', fontsize=22)
    ax.set_ylim(0, None)  # Set the lower limit of the y-axis to 0
    ax.yaxis.set_tick_params(labelsize=18)
    ax.xaxis.set_tick_params(labelsize=18)

    # Save the plot as a PNG file
    output_file = os.path.join(output_dir, f"{row['ID_contig']}_plot_phase.png")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close(fig)  # Close the figure after saving

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.csv output_directory")
    else:
        input_file = sys.argv[1]
        output_dir = sys.argv[2]

        # Create the output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        df = pd.read_csv(input_file, sep='\t')

        # Iterate over each row in the dataframe and create a plot for each
        for index, row in df.iterrows():
            plot_histogram(row, output_dir)
