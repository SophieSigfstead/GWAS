import os
import pandas as pd
import sys

def match_csv_tsv(directory_gwas_study_tsv, directory_1000_genomes, track_selection_list):
    chunksize = 100000
    output_path = './1000_genomes_csv_tsv_match'
    os.makedirs(output_path, exist_ok=True)
    original_result_file_path = os.path.join(output_path, 'combined_result_original.csv')
    sorted_result_file_path = os.path.join(output_path, 'combined_result_sorted.csv')

    # Load TSV data into memory
    tsv_data = pd.read_csv(directory_gwas_study_tsv, sep='\t')
    tsv_data['chromosome'] = pd.to_numeric(tsv_data['chromosome'], errors='coerce').astype('Int64')
    tsv_data['base_pair_location'] = pd.to_numeric(tsv_data['base_pair_location'], errors='coerce').astype('Int64')

    csv_columns_no_sad = ['alt', 'chr', 'pos', 'ref', 'snp']
    tsv_columns = ['beta', 'standard_error', 'effect_allele_frequency', 'p_value']

    # Columns to merge the files on
    csv_match_columns = ['chr', 'pos']
    tsv_match_columns = ['chromosome', 'base_pair_location']

    # Flag for headers being written or not
    headers_written = False

    with open(original_result_file_path, 'w') as f_out:
        csv_files = [f for f in os.listdir(directory_1000_genomes) if f.endswith('.csv') and f != 'targets.csv']
        total_files = len(csv_files)
        for i, csv_file in enumerate(csv_files):
            csv_file_path = os.path.join(directory_1000_genomes, csv_file)
            print(f"Processing file {i+1}/{total_files}: {csv_file}")
            for chunk_index, chunk in enumerate(pd.read_csv(csv_file_path, chunksize=chunksize)):
                print(f"Processing chunk {chunk_index+1} of {csv_file}")

                # Filter and calculate average SAD
                chunk_filtered = chunk[track_selection_list + csv_columns_no_sad].copy()
                chunk_filtered.loc[:, 'average_SAD'] = chunk_filtered[track_selection_list].mean(axis=1)  # Calculate average SAD

                merged_chunk = chunk_filtered.merge(
                    tsv_data,
                    left_on=csv_match_columns,
                    right_on=tsv_match_columns,
                    how='left'
                )

                # Remove duplicate columns from merged_chunk
                merged_chunk.drop(columns=['other_allele', 'chromosome', 'base_pair_location', 'effect_allele'], inplace=True)

                # Select and reorder the columns for output
                output_columns = track_selection_list + csv_columns_no_sad + tsv_columns + ['average_SAD']  # Include average_SAD
                merged_chunk = merged_chunk[output_columns]

                # Add header if required
                if not headers_written:
                    merged_chunk.to_csv(f_out, index=False, mode='a', header=True)
                    headers_written = True
                else:
                    merged_chunk.to_csv(f_out, index=False, mode='a', header=False)

    print(f"Original results written to {original_result_file_path}")

    # Sort by absolute value of average_SAD, and save as new file
    combined_data = pd.read_csv(original_result_file_path)
    combined_data['abs_average_SAD'] = combined_data['average_SAD'].abs()  # Calculate absolute value
    sorted_combined_data = combined_data.sort_values(by='abs_average_SAD')

    # Add ranking column
    sorted_combined_data['ranking'] = range(1, len(sorted_combined_data) + 1)

    sorted_combined_data.to_csv(sorted_result_file_path, index=False)
    print(f"Sorted results with ranking written to {sorted_result_file_path}")

    return original_result_file_path, sorted_result_file_path

def main(directory_gwas_study_tsv, directory_1000_genomes):
    original_path, sorted_path = match_csv_tsv(directory_gwas_study_tsv, directory_1000_genomes, ["SAD1", "SAD9"])
    return

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <path for GWAS tsv> <directory for 1000 genomes csvs>")
    else:
        directory_gwas_study_tsv = sys.argv[1]
        directory_1000_genomes = sys.argv[2]
        main(directory_gwas_study_tsv, directory_1000_genomes)
    sys.exit()
