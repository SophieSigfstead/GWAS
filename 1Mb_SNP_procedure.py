import pandas as pd
import sys

def main(path_to_combined_csv):

    # filter the CSV file to only include SNPs with non-null average_SAD
    df = pd.read_csv(path_to_combined_csv)
    df = df[df['average_SAD'].notna()]

    # create a csv file with the header 
    output_columns = ['SAD1', 'SAD9', 'alt', 'chr', 'pos', 'ref', 'snp', 'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'average_SAD', 'abs_average_SAD', 'ranking']
    result_df = pd.DataFrame(columns=output_columns)

    # set a threshold value, say THRESHOLD = 0.00002
    THRESHOLD = 0.00000

    # for the list of chromosomes in the tsv file (chr = 1 to 22), indicated by the column chr:
    for chr_num in range(1, 23):
        chr_df = df[df['chr'] == chr_num]
        chr_df = chr_df[(chr_df['average_SAD'] > THRESHOLD) | (chr_df['average_SAD'] < (-1*THRESHOLD))]
        chr_df = chr_df.sort_values(by='average_SAD', ascending=False).reset_index(drop=True)
        
        while not chr_df.empty:
            # get the top ranked SNP for that chromosome based on the average_SAD value, provided it is greater than THRESHOLD or less than -THRESHOLD. Add the full row to a csv file.
            top_snp = chr_df.iloc[0]
            result_df = pd.concat([result_df, pd.DataFrame([top_snp], columns=output_columns)], ignore_index=True)
            
            # let position_of_sig_snp = pos of the row just added to the csv.
            pos_sig_snp = top_snp['pos']
            
            # go down the list for the current chromosome only, and eliminate any rows that have (pos < position_of_sig_snp + 1000000) or (pos > position_of_sig_snp - 1000000)
            chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + 100000) & (chr_df['pos'] > pos_sig_snp - 100000))]
            chr_df = chr_df.reset_index(drop=True)

    # save the csv file of the rows found in the above loop
    result_df.to_csv('./filtered_snps.csv', index=False)
    
    # print out the length of the saved csv file
    print(f"Length of the saved CSV file: {len(result_df)}")
    
    # using another file located in ./sorted_GWAS_summary_data/GCST90277450_sorted_significant.csv
    other_file_path = './sorted_GWAS_summary_data/GCST90277450_sorted_significant.csv'
    other_df = pd.read_csv(other_file_path, sep='\t')
    
    # print out how many rows share the same value for chr (or chromosome in the GCST90277450_sorted_significant.csv)
    common_chr_count = len(set(result_df['chr']).intersection(set(other_df['chromosome'])))
    print(f"Number of common chromosomes: {common_chr_count}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path for combined csv>")
    else:
        path_to_combined_csv = sys.argv[1]
        main(path_to_combined_csv)
    sys.exit()
