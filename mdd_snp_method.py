import sys
import pandas as pd 

def main(path_to_tsv):

    # given the path to the csv
    df = pd.read_csv(path_to_tsv, sep='\t')

    # create a csv file with the header 
    # chromosome	base_pair_location	effect_allele	other_allele	beta	standard_error	effect_allele_frequency	p_value	ranking_of_p_value
    output_columns = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'ranking_of_p_value']
    
    # this will store the result
    result_df = pd.DataFrame(columns=output_columns)
    
    for chr_num in range(1, 23):

        chr_df = df[df['chromosome'] == chr_num]
        chr_df = chr_df.sort_values(by='p_value', ascending=True).reset_index(drop=True)
        
        while not chr_df.empty:
            # get the top ranked SNP for that chromosome based on pvalue, provided it is greater than THRESHOLD or less than -THRESHOLD. Add the full row to a csv file.
            top_snp = chr_df.iloc[0]
            result_df = pd.concat([result_df, pd.DataFrame([top_snp], columns=output_columns)], ignore_index=True)
            
            # let position_of_sig_snp = pos of the row just added to the csv.
            pos_sig_snp = top_snp['base_pair_location']
            
            # go down the list for the current chromosome only, and eliminate any rows that have (pos < position_of_sig_snp + 1000000) or (pos > position_of_sig_snp - 1000000)
            chr_df = chr_df[~((chr_df['base_pair_location'] < pos_sig_snp + 1000000) & (chr_df['base_pair_location'] > pos_sig_snp - 1000000))]
            chr_df = chr_df.reset_index(drop=True)
    
    # save the csv file of the rows found in the above loop
    result_df.to_csv('./mdd_sig_snps.csv', index=False)
    
    # print out the length of the saved csv file
    print(f"Length of the saved mdd_sig_snps file: {len(result_df)}")

    return

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path for tsv>")
    else:
        path_to_tsv = sys.argv[1]
        main(path_to_tsv)
    sys.exit()