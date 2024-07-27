import pandas as pd
import os
import sys
from intervaltree import IntervalTree

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

    exon_regions = pd.read_csv('./genome_assembly/exon_regions.csv')

    # take out single location exons from the intervals, so can build interval tree
    exon_regions['diff'] = exon_regions['end'] - exon_regions['start']
    single_locations =  exon_regions[exon_regions['diff']==0]
    single_locations_list = list(single_locations['start'])

    # define interval tree
    exon_regions = exon_regions[exon_regions['diff']>=1]
    interval_tree = IntervalTree()
    for start, end in zip(exon_regions['start'], exon_regions['end']):
        interval_tree[start:end+1] = True
        
    def is_in_intervals(pos):
        in_intervals = bool(interval_tree[pos])
        in_single_regions = pos in single_locations_list
        return in_intervals or in_single_regions

    # Read the first file to initialize the in_coding_region column
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))

    # Ensure 'pos' is aligned for the interval check
    first_file_df['in_coding_region'] = first_file_df['pos'].apply(lambda pos: is_in_intervals(pos))

    print(len(first_file_df))
    
    coding_region_set = set(first_file_df[first_file_df['in_coding_region']]['snp'])
    print('length of coding region set')
    print(len(coding_region_set))

    for file in combined_files:
        print("processing:", file)
        df = pd.read_csv(os.path.join(directory_gwas_combined_files, file), low_memory = False)
        df = df[df['snp'].notna() & ~df['snp'].isin(coding_region_set)]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]  
        upper_threshold, lower_threshold = threshold_calculator(df, sd, track_col)
   
        above_threshold = df[(df[track_col] > upper_threshold) | (df[track_col] < lower_threshold)]
        below_threshold = df[(df[track_col] <= upper_threshold) & (df[track_col] >= lower_threshold)]
        
        threshold_set.update(above_threshold['snp'])
        non_threshold_set.update(below_threshold['snp'])
    print("sd:", sd)
    
    print('n snps total', len(first_file_df))

    threshold_df = first_file_df[first_file_df['snp'].isin(threshold_set)]

    non_threshold_df = first_file_df[first_file_df['snp'].isin(non_threshold_set)]
    
    common_marker_names = set(threshold_df['snp']).intersection(set(non_threshold_df['snp']))

    if common_marker_names:

        non_threshold_df = non_threshold_df[~non_threshold_df['snp'].isin(common_marker_names)]
    
    coding_non_coding_crossover = set(threshold_set.intersection(coding_region_set)).union(non_threshold_set.intersection(coding_region_set))
    
    if coding_non_coding_crossover:
        print("CODING NON CODING crossover")


    print('n snps in sets')
    print(len(threshold_df) + len(non_threshold_df) + len(coding_region_set))

    
    # Modification: the n_snps_prev should be equal to the length of the threshold and non-threshold df together
    n_snps_prev = len(threshold_df) + len(non_threshold_df) + len(coding_region_set)
    print("n snps prev", n_snps_prev)

    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region'])]
    
    n_snps_filtered = len(threshold_df) + len(coding_region_df)
    print("n snps filtered", n_snps_filtered)
    
    new_pvalue = 5e-8 * (n_snps_prev / n_snps_filtered)
    print("New p-value threshold: ", new_pvalue)

    threshold_df = threshold_df[(threshold_df['p_value'].notna()) & (threshold_df['p_value'] < new_pvalue)]
    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region']) & (first_file_df['p_value'] < new_pvalue)]
    
    threshold_and_coding = pd.concat([threshold_df, coding_region_df], ignore_index=True)

    result_df = pd.DataFrame()
    
    for chr_num in range(1, 23):
        chr_df = threshold_and_coding[threshold_and_coding['chr'] == chr_num ].sort_values(by='p_value').reset_index(drop=True)
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = pd.concat([result_df, pd.DataFrame([top_snp])], ignore_index=True)
            pos_sig_snp = top_snp['pos']
            chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + 1000000) & (chr_df['pos'] > pos_sig_snp - 1000000))]
            chr_df = chr_df.reset_index(drop=True)
    
    
    print('Sig snps coding:', len(coding_region_df))

    os.makedirs('./filtered_snps_gwas_1_v2', exist_ok=True) 
    result_df.to_csv(f'./filtered_snps_gwas_1_v2/filtered_snps_gwas_1_threshold={sd}.csv', index=False)


    print('number of snps found', len(result_df))
    
    # GWAS 2
    # compare SNPs with the significant SNPs file
    #sig_snps_df = pd.read_csv(directory_gwas_sig_snps)
    #common_snps = set(result_df['snp']).intersection(set(sig_snps_df['Marker_Name']))
    #print(f'Number of common SNPs !!: {len(common_snps)}')

    #GWAS #1
    mdd_sig_snps = pd.read_csv('mdd_sig_snps.csv')
    matching_snps = pd.merge(result_df, mdd_sig_snps, left_on=['chr', 'pos'], right_on=['chromosome', 'base_pair_location'])
    print(f"Number of matching rows with mdd_sig_snps.csv: {len(matching_snps)}")

    result_df['right_border'] = result_df['pos'] + 1000000
    result_df['left_border'] = result_df['pos'] - 1000000

    result_df.to_csv(f'./filtered_snps_gwas_1_v2/filtered_snps_test.csv', index=False)

    mdd_sig_snps['right_border'] = mdd_sig_snps['base_pair_location'] + 1000000
    mdd_sig_snps['left_border'] = mdd_sig_snps['base_pair_location'] - 1000000

    overlap_count = 0
    print('YEP')
    # Iterate through each interval in df1
    for i, row1 in mdd_sig_snps.iterrows():
        left1, right1 = row1['left_border'], row1['right_border']
        chr1  = row1['chromosome']
        
        # Check for overlap with each interval in df2
        for j, row2 in result_df.iterrows():
            left2, right2 = row2['left_border'], row2['right_border']
            chr2  = row2['chr']
            
            # Check if the intervals overlap
            if (left1 <= right2 and right1 >= left2) and (chr1 == chr2):
                overlap_count += 1
                break  # Exit the loop once an overlap is found for the current interval in df1
    
                
    print("overlapping loci: " ,overlap_count)

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
