import sys
import pandas as pd

# Step 1: Process the contig file and add RSState
def process_contig_file(input_file):
    df_contig = pd.read_csv(input_file, sep='\t')
    df_contig['RSState'] = df_contig['p_value'].apply(lambda x: 'RS+P+' if x < 0.05 else 'RS+P-')
    return df_contig

# Step 2: Process the second file and append it to the contig output
def process_and_merge_files(contig_df, second_file, output_file):
    df_second = pd.read_csv(second_file, sep='\t')
    second_output_df = pd.DataFrame()
    second_output_df['contig'] = df_second['tag']
    second_output_df['ID_contig'] = df_second['ID']
    for col in contig_df.columns[2:-1]:  # Skip the 'contig', 'ID_contig', and 'RSState' columns
        second_output_df[col] = 'NA'
    second_output_df['RSState'] = 'RS-'
    combined_df = pd.concat([contig_df, second_output_df], ignore_index=True)
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"File saved successfully as {output_file}")

# Main function
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <contig_file> <second_file> <output_file>")
        sys.exit(1)
    contig_file = sys.argv[1]
    second_file = sys.argv[2]
    output_file = sys.argv[3]

    # Step 1: Process the contig file
    contig_df = process_contig_file(contig_file)

    # Step 2: Process the second file and merge with the contig data
    process_and_merge_files(contig_df, second_file, output_file)
