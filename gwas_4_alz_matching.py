import os
import pandas as pd
import sys
import ast
import shutil
import gc

def clear_directory(directory_path):
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    os.makedirs(directory_path)

def match_csv_txt(gwas_directory, directory_1000_genomes, track_list):
   
    chunksize = 100000
    output_path = './gwas_4_alz_matching'
    clear_directory(output_path)
    gc.collect()
    text_data = pd.DataFrame()
    # Load text file data into memory and filter relevant columns
    text_data = pd.read_csv(gwas_directory, sep='\t' ,names = ['chr',	'PosGRCh37',	'testedAllele',	'otherAllele',	'z',	'p',	'N'])
    text_data = text_data[['chr',	'PosGRCh37'	,'testedAllele',	'otherAllele',	'z',	'p',	'N']]
    print(len(text_data))
    sum_stats_2019 = pd.read_csv('gwas_4_alz_sum_stats_2019_Kunkle', sep = '\t')
    text_with_snp_ids = pd.merge(text_data, pd.DataFrame(sum_stats_2019[['SNP', 'CHR', "BP"]]), left_on =['chr','PosGRCh37' ], right_on=['CHR', "BP"], how='left')
    print(len(text_with_snp_ids))

    #significant_with_snp_ids that are nan look up in another dataframe called gwas_1_matching
    gwas_1_matching = pd.read_csv('gwas_1_matching/result_SAD0.csv')
    nan_snp_ids = text_with_snp_ids[text_with_snp_ids['SNP'].isna()] 
    nan_snp_ids = pd.merge(nan_snp_ids, gwas_1_matching[['chr', 'pos', 'snp']], left_on=['chr', 'PosGRCh37'], right_on=['chr', 'pos'], how='left')
    text_with_snp_ids.loc[text_with_snp_ids['SNP'].isna(), 'SNP'] = nan_snp_ids['snp']
    text_data = text_with_snp_ids[text_with_snp_ids['SNP'].notna()]
   
    text_data['SNP'] = text_data['SNP'].astype(str)
    print(len(text_data))

    csv_columns_no_sad = ['alt', 'chr_x', 'pos', 'ref', 'snp']
    csv_match_columns = ['snp']

    # Dictionary to keep track of which files have headers written
    headers_written = {sad: False for sad in track_list}

    # Process each CSV file in the 1000 genomes directory
    csv_files = [f for f in os.listdir(directory_1000_genomes) if f.endswith('.csv') and f != 'targets.csv']
    total_files = len(csv_files)
    for i, csv_file in enumerate(csv_files):
        csv_file_path = os.path.join(directory_1000_genomes, csv_file)
        print(f"Processing file {i+1}/{total_files}: {csv_file}")
        for chunk_index, chunk in enumerate(pd.read_csv(csv_file_path, chunksize=chunksize)):
            print(f"Processing chunk {chunk_index+1} of {csv_file}")
            chunk['snp'] = chunk['snp'].astype(str)
           
            if 'snp' not in chunk.columns:
                print(f"'snp' column not found in chunk {chunk_index+1} of {csv_file}")
                continue

            # Merge the entire chunk once with the text data
            merged_chunk = chunk.merge(text_data, how='left', left_on=['snp'], right_on=['SNP'])
            #print(merged_chunk[['snp', 'SNP']].head())

            for track in track_list:
                sad_column = f'SAD{track}'
                if sad_column not in merged_chunk.columns:
                    print(f"{sad_column} not found in chunk {chunk_index+1} of {csv_file}")
                    continue
                
                chunk_filtered_sad = merged_chunk[[sad_column] + csv_columns_no_sad + ['chr_y',	'PosGRCh37'	,'testedAllele',	'otherAllele',	'z',	'p',	'N', 'SNP']].reset_index(drop=True)

                # Define the output file path for the current track
                track_output_file_path = os.path.join(output_path, f'result_SAD{track}.csv')

                # Write the merged chunk to the appropriate file
                if not headers_written[track]:
                    chunk_filtered_sad.to_csv(track_output_file_path, index=False, mode='a', header=True)
                    headers_written[track] = True
                else:
                    chunk_filtered_sad.to_csv(track_output_file_path, index=False, mode='a', header=False)
    
        # Load the first result file to get the matched MarkerNames
    new_df = pd.read_csv(os.path.join(output_path, 'result_SAD0.csv'))
    matched_markernames = set(new_df['SNP'])
    #print(matched_markernames)

        # Remove matched MarkerNames from text_data
    remaining_text_data = text_data[~text_data['SNP'].isin(matched_markernames)].reset_index()
    print(len(remaining_text_data))

    
        # Add remaining text rows to each combined file with null values for the columns from the CSV files
    output_path = './gwas_4_alz_matching'
    for track in track_list:
        track_output_file_path = os.path.join(output_path, f'result_SAD{track}.csv')
        remaining_text_data_with_nulls = remaining_text_data.copy()
        for col in [f'SAD{track}'] + csv_columns_no_sad:
            remaining_text_data_with_nulls[col] = pd.NA    
        if not remaining_text_data_with_nulls.empty:
            remaining_text_data_with_nulls = remaining_text_data_with_nulls[[f'SAD{track}'] + csv_columns_no_sad + ['chr',	'PosGRCh37'	,'testedAllele',	'otherAllele',	'z',	'p',	'N', 'SNP']]
            remaining_text_data_with_nulls.to_csv(track_output_file_path, index=False, mode='a', header=False)

    print("Processing complete.")

def main(directory_gwas_study_text, directory_1000_genomes, track_list):
    track_list = ast.literal_eval(track_list)  # Convert string representation of list to actual list
    match_csv_txt(directory_gwas_study_text, directory_1000_genomes, track_list)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <path for GWAS text> <directory for 1000 genomes csvs> <track list e.g. [1,2,3]>")
    else:
        directory_gwas_study_text = sys.argv[1]
        directory_1000_genomes = sys.argv[2]
        track_list = sys.argv[3]
        main(directory_gwas_study_text, directory_1000_genomes, track_list)
    sys.exit()