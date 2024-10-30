import sys
import pandas as pd
import os 
from intervaltree import IntervalTree
import shutil

def threshold_calculator(df, sd, track):
    std_dev = df[track].std()
    mean = df[track].mean()

    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)

    return upper_threshold, lower_threshold

def clear_directory(directory_path):
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    os.makedirs(directory_path)

def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, key):

    sd = float(sd)

    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]

    ## STEP 1: Define the coding SNPs list
    exon_regions = pd.read_csv('../genome_assembly/exon_regions_v2.csv')
    # Read the first file to initialize the in_coding_region column
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))
    '''
    # create interval trees for each exon in exon regions 
    trees = {}
    for _, row in exon_regions.iterrows():
        chr = row['chr']
        if chr not in trees:
            trees[chr] = IntervalTree()
            print(f'chromosome {chr} has been added')
        else:
            trees[chr][row['start']:row['end']+1] = True
    
    # check if a SNP is within any exon region
    def is_within_exon(row):
        chrom = row['chr']
        pos = row['pos']
        if bool(trees[chrom][pos]):
            return True
        return False
    
    first_file_df = first_file_df[first_file_df['snp'].notna()]
    # Apply the function to check for overlaps
    first_file_df['in_coding_region'] = first_file_df.apply(is_within_exon, axis=1)
    
    coding_region_set = list(set(first_file_df[first_file_df['in_coding_region'] & first_file_df['p_value'].notna()]['snp']))
    # Convert to DataFrame
    coding_region_df = pd.DataFrame(coding_region_set, columns=['snp'])

    # Save as CSV
    coding_region_df.to_csv('./coding_region_set.csv', index=False)

    '''
    coding_region_df = pd.read_csv('./coding_region_set.csv')
    first_file_df['in_coding_region'] = first_file_df['snp'].isin(list(coding_region_df['snp']))
    coding_region_set = list(set(first_file_df[first_file_df['in_coding_region'] & first_file_df['p_value'].notna()]['snp']))
    

    # Create a directory to store the results of each analysis 
    clear_directory(f"./gwas_1_leading_SNPs_by_track_random_results/gwas_1_leading_SNPs_by_track_random_{key}/sd={sd}")

    ## STEP 2: For each file in the combined files directory, perform the GWAS analysis
    for file in combined_files:
        print(f"Starting {file}")
        threshold_set = set()
        non_threshold_set = set()
        
        df = pd.read_csv(os.path.join(directory_gwas_combined_files, file), low_memory = False)
        df = df[df['snp'].notna() & df['p_value'].notna()].reset_index()
    
        coding_region_df = df[df['snp'].isin(coding_region_set)]

        filt_df = df[~(df['snp'].isin(coding_region_set))]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]  
        
        upper_threshold, lower_threshold = threshold_calculator(df, sd, track_col)
        above_threshold = filt_df[(filt_df[track_col] > upper_threshold) | (filt_df[track_col] < lower_threshold)]
        below_threshold = filt_df[(filt_df[track_col] <= upper_threshold) & (filt_df[track_col] >= lower_threshold)]
        
        threshold_set.update(above_threshold['snp'])
        non_threshold_set.update(below_threshold['snp'])

        threshold_df = filt_df[filt_df['snp'].isin(threshold_set)]
        non_threshold_df = filt_df[filt_df['snp'].isin(non_threshold_set)]
        
        common_marker_names = threshold_set.intersection(non_threshold_set)
        if common_marker_names:
            print("common marker!")
            non_threshold_df = non_threshold_df[~non_threshold_df['snp'].isin(common_marker_names)]
        
        coding_and_threshold_df = pd.concat([threshold_df, coding_region_df], ignore_index=True)

        n_prev = len(df)
        print(f"n_prev: {n_prev}")
        n_filtered = len(coding_and_threshold_df)
        print(f"n_filtered: {n_filtered}")
        print(f" len non-threshold df", len(non_threshold_df))
        print("len threshold and coding together == n_filtered"  , len(threshold_df)+len(coding_region_df))

        new_pvalue = 5e-8 * ( n_prev / n_filtered )

        coding_and_threshold_df = coding_and_threshold_df[coding_and_threshold_df['p_value']< new_pvalue]    
    
        result_df = pd.DataFrame()

        for chr_num in range(1, 23):
            chr_df = coding_and_threshold_df[coding_and_threshold_df['chr'] == chr_num ].sort_values(by='p_value').reset_index(drop=True)
            while not chr_df.empty:
                top_snp = chr_df.iloc[0]
                result_df = pd.concat([result_df, pd.DataFrame([top_snp])], ignore_index=True)
                pos_sig_snp = top_snp['pos']
                chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + 1000000) & (chr_df['pos'] > pos_sig_snp - 1000000))]
                chr_df = chr_df.reset_index(drop=True)   
    
        # store the new GWAS results
        result_df['left_border'] = result_df['pos'] - 1000000
        result_df['right_border'] = result_df['pos'] + 1000000
        print(result_df.head())
        result_df.to_csv(f'./gwas_1_leading_SNPs_by_track_random_results/gwas_1_leading_SNPs_by_track_random_{key}/sd={sd}/filtered_snps_gwas_1_track={track_col}_sd={sd}.csv', index=False)

        # compute the comparison with the original study, and then store these results in a dataframe
        overlap_df = pd.DataFrame(columns =['track_col', 'matching_snps', "matching_loci", "new_coding_snps"])

        # compute intersection of snp list
        mdd_sig_snps = pd.read_csv('../gwas_1_intermediate_files/mdd_sig_snps.csv')
        
        matching_snps = pd.merge(result_df, mdd_sig_snps, left_on=['chr', 'pos'], right_on=['chromosome', 'base_pair_location'])

        # compute number of overlapping loci 
        mdd_sig_snps['right_border'] = mdd_sig_snps['base_pair_location'] + 1000000
        mdd_sig_snps['left_border'] = mdd_sig_snps['base_pair_location'] - 1000000

        overlap_count = 0
        
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

        # compute number of new coding snps
        coding_snps_result_df =  result_df[result_df['snp'].isin(coding_region_set)]
    
        # store the results
        overlap_df = overlap_df._append({'track_col': track_col, 'matching_snps': len(matching_snps), "matching_loci": overlap_count, "new_coding_snps": len(coding_snps_result_df)}, ignore_index=True)

        overlap_df.to_csv(f'./gwas_1_leading_SNPs_by_track_random_results/gwas_1_leading_SNPs_by_track_random_{key}/sd={sd}/filtered_snps_gwas_1_track={track_col}_sd={sd}_overlap_df.csv', index=False)

        print(f"Analysis for {track_col} complete.")

    return 


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        directory_gwas_sig_snps = sys.argv[2]
        sd = sys.argv[3]
        key = sys.argv[4]
        main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, key)