import pandas as pd
import os
import sys

def threshold_calculator(df, sd, track):
    df_with_sad = df[(df[track].notna()) & (df[track].abs() < 1) & (df[track].abs() >= 0)]
  
    std_dev = df_with_sad[track].std()
    mean = df_with_sad[track].mean()

    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)

    return upper_threshold, lower_threshold

def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd):
    sd = float(sd)
    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]
    
    threshold_set = set()
    non_threshold_set = set()
    
    for file in combined_files:
        print("processing:", file)
        df = pd.read_csv(os.path.join(directory_gwas_combined_files, file))
        df = df[df['snp'].notna()]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]  # Assumes only one SAD column
        upper_threshold, lower_threshold = threshold_calculator(df, sd, track_col)
        
        above_threshold = df[(df[track_col] > upper_threshold) | (df[track_col] < lower_threshold)]
        below_threshold = df[(df[track_col] <= upper_threshold) & (df[track_col] >= lower_threshold)]
        
        threshold_set.update(above_threshold['snp'])
        non_threshold_set.update(below_threshold['snp'])
    
    print("before drop", len(threshold_set))

    # Read the first file to get SNP data
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))

    threshold_df = first_file_df[first_file_df['snp'].isin(threshold_set)]

    non_threshold_df = first_file_df[first_file_df['snp'].isin(non_threshold_set)]
    
    common_marker_names = set(threshold_df['snp']).intersection(set(non_threshold_df['snp']))

    if common_marker_names:
        # Remove common markers from non_threshold_df
        non_threshold_df = non_threshold_df[~non_threshold_df['snp'].isin(common_marker_names)]
    
    n_snps_prev = len(first_file_df[first_file_df.filter(like='SAD').notna().any(axis=1)])
    print("n snps prev", n_snps_prev)

    n_snps_filtered = len(threshold_df)
    print("n snps filtered", n_snps_filtered)
    
    new_pvalue = 5e-8 * (n_snps_prev / n_snps_filtered)
    print("New p-value threshold: ", new_pvalue)
    
    threshold_df = threshold_df[(threshold_df['P'].notna()) & (threshold_df['P']< new_pvalue)]
    
    result_df = pd.DataFrame()
    
    for chr_num in range(1, 23):
        chr_df = threshold_df[threshold_df['chr'] == chr_num].sort_values(by='P').reset_index(drop=True)
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = pd.concat([result_df, pd.DataFrame([top_snp])], ignore_index=True)
            pos_sig_snp = top_snp['pos']
            chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + 1000000) & (chr_df['pos'] > pos_sig_snp - 1000000))]
            chr_df = chr_df.reset_index(drop=True)
    
    os.makedirs('./filtered_snps_gwas_2', exist_ok=True)
    
    result_df.to_csv(f'./filtered_snps_gwas_2/filtered_snps_gwas_2_threshold={sd}.csv', index=False)

    # Compare SNPs with the significant SNPs file
    sig_snps_df = pd.read_csv(directory_gwas_sig_snps)
    common_snps = set(result_df['snp']).intersection(set(sig_snps_df['Marker_Name']))
    print(f'Number of common SNPs: {len(common_snps)}')
    return
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        directory_gwas_sig_snps = sys.argv[2]
        sd = sys.argv[3]
        main(directory_gwas_combined_files, directory_gwas_sig_snps, sd)
    sys.exit()
