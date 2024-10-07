import pandas as pd
import os
import sys
from intervaltree import IntervalTree
import subprocess

def threshold_calculator(df, sd, track):
    df_with_sad = df[(df[track].notna()) & (df[track].abs() < 1) & (df[track].abs() >= 0)]
  
    std_dev = df_with_sad[track].std()
    mean = df_with_sad[track].mean()

    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)

    return upper_threshold, lower_threshold

def merge_loci(df, distance_threshold=250000):
    merged = []
    df = df.sort_values(by=['chr', 'left_border']).reset_index(drop=True)  # Sort by chromosome and left border
    
    while not df.empty:
        current = df.iloc[0]  # Take the first row
        df = df.iloc[1:]  # Remove the first row from the dataframe
        
        # Find overlapping or nearby rows within the threshold
        overlapping_rows = df[
            (df['chr'] == current['chr']) &
            (df['left_border'] <= current['right_border'] + distance_threshold) &
            (df['right_border'] >= current['left_border'] - distance_threshold)
        ]
        
        # If there are overlaps, merge them
        if not overlapping_rows.empty:
            # Update the current interval
            print("overlaping rows", len(overlapping_rows))
            overlapping_rows = overlapping_rows.sort_values(['p'], ascending = True)

            min_left_border = min(current['left_border'], overlapping_rows['left_border'].min())
            max_right_border = max(current['right_border'], overlapping_rows['right_border'].max())
            
            # Get the row with the lowest p-value
            min_p_row = overlapping_rows.loc[0]
            
            # Create the merged interval
            merged_interval = {
                'pos': min_p_row['pos'],
                'chr': current['chr'],
                'p': min_p_row['p'],
                'snp': min_p_row['snp'],
                'left_border': min_left_border,
                'right_border': max_right_border
            }
            
            # Remove the overlapping rows from the dataframe
            df = df.drop(overlapping_rows.index).reset_index(drop=True)
        else:
            # If no overlap, the current interval is finalized
            merged_interval = current.to_dict()
        
        merged.append(merged_interval)
    
    return pd.DataFrame(merged)

def main(directory_gwas_combined_files, sd):

    sd = float(sd) 

    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]
    
    threshold_set = set()
   
    non_threshold_set = set()

    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]), low_memory = False)
    print("first_file_df", first_file_df.head())

    non_match_list = list(first_file_df[(first_file_df['snp'].isna()) & (first_file_df['p'].notna())]['snp'])
    print("number of snps", len(first_file_df))
    print("number of snps in gwas not in 1KG: ", len(non_match_list))
    non_gwas_list = list(first_file_df[(first_file_df['p'].isna()) & (first_file_df['snp'].notna())]['snp'])
    print("number of snps in 1kg not in gwas: ", len(non_gwas_list))
    
    '''
    exon_regions = pd.read_csv('./genome_assembly/exon_regions_v2.csv')
    
    # create interval trees for each exon in exon regions 
    trees = {}
    for _, row in exon_regions.iterrows():
        chr = row['chr']
        if chr not in trees:
            #print(f"{chr} was added")
            trees[chr] = IntervalTree()
        trees[chr][row['start']:row['end']+1] = True

    # check if a SNP is within any exon region
    def is_within_exon(row):
        chrom = row['chr_x']
        pos = row['pos']
        if chrom in trees:
            return bool(trees[chrom][pos])
        return False

    # Apply the function to check for overlaps
    first_file_df['in_coding_region'] = first_file_df.apply(is_within_exon, axis=1)
   
    coding_region_set = set(first_file_df[first_file_df['in_coding_region']]['snp'])
    first_file_df[first_file_df['in_coding_region']].to_csv("gwas_4_alz_intermediate_files/coding_snps.csv") 
    '''
    coding_snps_df = pd.read_csv("gwas_4_alz_intermediate_files/coding_snps.csv")
    coding_region_set = set(coding_snps_df['snp'])
    print('length of coding region set')
    print(len(coding_region_set))

    for file in combined_files:
        print("processing:", file)
        df = pd.read_csv(os.path.join(directory_gwas_combined_files, file), low_memory = False)
        df = df[df['snp'].notna()]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]  
        upper_threshold, lower_threshold = threshold_calculator(df, sd, track_col)
   
        above_threshold = df[(df[track_col] > upper_threshold) | (df[track_col] < lower_threshold)]
        below_threshold = df[(df[track_col] <= upper_threshold) & (df[track_col] >= lower_threshold)]
        
        threshold_set.update(above_threshold['snp'])
        non_threshold_set.update(below_threshold['snp'])
    print("sd:", sd)
    
    print('n snps total', len(first_file_df))

    threshold_df = first_file_df[(first_file_df['snp'].isin(threshold_set)) & (first_file_df['snp'].notna()) & (first_file_df['p'].notna()) & (~first_file_df['snp'].isin(coding_region_set))]

    non_threshold_df = first_file_df[first_file_df['snp'].isin(non_threshold_set) &  (first_file_df['snp'].notna()) & (first_file_df['p'].notna()) & ~first_file_df['snp'].isin(coding_region_set)]
    
    # take out snps that met another threshold from non-threshold-df
    common_marker_names = set(threshold_df['snp']).intersection(set(non_threshold_df['snp']))
    if common_marker_names:
        non_threshold_df = non_threshold_df[~non_threshold_df['snp'].isin(common_marker_names)]

    non_threshold_df.to_csv(f'gwas_4_alz_intermediate_files/gwas_4_non_threshold_snps_sd={sd}.csv') 
    coding_region_df = first_file_df[(first_file_df['snp'].isin(coding_region_set)) & (first_file_df['p'].notna())]
    # Modification: the n_snps_prev should be equal to the length of the threshold and non-threshold df together
    n_snps_prev = len(threshold_df) + len(non_threshold_df) + len(coding_region_df)
    print("n snps prev - should equal number in the study: ", n_snps_prev)
    
    n_snps_filtered = len(threshold_df) + len(coding_region_df) 
    print("n snps after filtering: ", n_snps_filtered)
    
    new_pvalue = 5e-8 * (n_snps_prev / n_snps_filtered)

    print("New p-value threshold: ", new_pvalue)
    
    threshold_df['p'] = pd.to_numeric(threshold_df['p'], errors='coerce')
    print(threshold_df[threshold_df['p'].isna()])

    coding_region_df['p'] = pd.to_numeric(coding_region_df['p'],  errors='coerce')
    print(coding_region_df[coding_region_df['p'].isna()])
    #coding_region_df = coding_region_df[coding_region_df['P'] < new_pvalue]
    
    threshold_and_coding = pd.concat([threshold_df, coding_region_df], ignore_index=True)

    threshold_and_coding_sorted = threshold_and_coding.sort_values('p').reset_index()
    threshold_and_coding_sorted = threshold_and_coding_sorted[threshold_and_coding_sorted['snp'].notna()]
    threshold_and_coding_sorted.to_csv(f'gwas_4_alz_intermediate_files/gwas_4_threshold_and_coding_snps_sd={sd}.csv') 
    
    
    threshold_and_coding_sorted = pd.read_csv(f'gwas_4_alz_intermediate_files/gwas_4_threshold_and_coding_snps_sd={sd}.csv')
    print("threshold_and_coding_sorted")
    print(threshold_and_coding_sorted.head())
    threshold_and_coding_sorted = threshold_and_coding_sorted[threshold_and_coding_sorted['p'] < new_pvalue]
    
    
    plink_path = "../tools/plink2/plink2"
    test_file_path = './1KG_genotypes/all_phase3_bin_eur_no_duplicates'

    all_variants = pd.DataFrame(threshold_and_coding_sorted['SNP']).rename(columns={'SNP': 'ID'})

    all_variants.to_csv('./gwas_4_alz_intermediate_files/all_variants_list', sep = "\t", index=False)
    result_df = pd.DataFrame(columns = ['pos', 'chr', 'p', 'snp', 'left_border', 'right_border'])

    
    for _, row in threshold_and_coding_sorted.iterrows():
        pos = int(row['pos'])
        chr = int(row['chr_x'])
        p = row["p"]
        snp = row["SNP"]
        start = int(pos - 1500000)
        end = int(pos + 1500000)


        plink_make_bed_command = [
        f'{plink_path}', 
        '--bfile', test_file_path,  
        '--chr', str(chr),  
        '--from-bp', str(start),  
        '--to-bp', str(end),  
        '--make-bed',  
        '--out', f'./gwas_4_alz_intermediate_files/chr_bed_files/chr_bed_file_chr={chr}_pos={pos}'
        ]
        subprocess.run(plink_make_bed_command, check=True)

        plink_ld_command = [
            f'{plink_path}', 
            "--bfile", f'./gwas_4_alz_intermediate_files/chr_bed_files/chr_bed_file_chr={chr}_pos={pos}',  
            "--ld-window", "99999",
            "--ld-window-kb", "3000",
            "--ld-window-r2", "0.6",
            "--r2-unphased",
            "--ld-snp-list",  f'{'./gwas_4_alz_intermediate_files/all_variants_list'}',
            "--out", f"./gwas_4_alz_intermediate_files/ld_output_files/ld_loci_output_chr={chr}_pos={pos}"
        ]
        subprocess.run(plink_ld_command)

        def define_locus(snp, chr, pos):
            ld_data = pd.read_csv(f'./gwas_4_alz_intermediate_files/ld_output_files/ld_loci_output_chr={chr}_pos={pos}.vcor', sep = '\t')
            correlated_snps = ld_data[(ld_data['ID_A'] == snp) & (ld_data['UNPHASED_R2'] > 0.6)]
            positions = correlated_snps['POS_B'].tolist()
            if ld_data[(ld_data['ID_A'] == snp)].empty:
                print("snp doesn't exist in ld data")
            if positions:
                locus_start = min(positions)
                locus_end = max(positions)
            else:
                locus_start = locus_end = pos
            return locus_start, locus_end

        # define locus based on ld information
        left_border, right_border = define_locus(snp, chr, pos)
        result_df=result_df._append({"chr":chr, "pos": pos, "snp": snp, "p": p, "left_border": left_border, "right_border": right_border}, ignore_index=True)
        print(result_df.head())
        print(len(result_df))

        result_df = result_df.sort_values(by=['chr', 'pos'])
        result_df.to_csv(f"./gwas_4_alz_intermediate_files/result_df_sd={sd}.csv")
        
        merged_result_df = merge_loci(result_df)
        merged_result_df = merge_loci(merged_result_df)
        
        merged_result_df.to_csv(f'./gwas_4_alz_result_files/filtered_snps_sd={sd}.csv')
        print("length of merged df " , len(merged_result_df))
    
    ground_truth_snps = pd.read_csv('./gwas_4_alz_result_files/filtered_snps_sd=0.0.csv')
    clumped_results_with_merged_borders = merged_result_df.copy()

    def count_overlaps(df1, df2):
        
        overlap_count = 0
    
        # Iterate over all pairs of intervals
        for _, row1 in df1.iterrows():
            chr1 = int(row1['chr'])
            for _, row2 in df2.iterrows():
                chr2 = int(row2['chr'])
                # Check if the intervals overlap on the same chromosome
                if chr1==chr2 and row1['right_border'] >= row2['left_border'] and row1['left_border'] <= row2['right_border']:
                    overlap_count += 1
                    break; 
        return overlap_count

    print("number of snps found with our method", len(clumped_results_with_merged_borders))

    filtered_in_coding_region = clumped_results_with_merged_borders[clumped_results_with_merged_borders['snp'].isin(coding_region_set)]

    print("number of coding snps" , len(filtered_in_coding_region))
    print(filtered_in_coding_region)
    filtered_in_coding_region.to_csv(f'gwas_4_alz_intermediate_files/our_coding_snps_sd={sd}')

    ground_truth_coding = ground_truth_snps[ground_truth_snps['snp'].isin(coding_region_set)]
    print("number of coding snps originally in ground truth: ", len(ground_truth_coding))
    print(ground_truth_coding)
    ground_truth_coding.to_csv(f'gwas_4_alz_intermediate_files/ground_truth_coding_snps_sd={sd}')

    overlap_count = count_overlaps(clumped_results_with_merged_borders, ground_truth_snps)

    print(f"Number of overlapping loci: {overlap_count}")

    common_snps = set(ground_truth_snps['snp']).intersection(clumped_results_with_merged_borders['snp'])

    print(f"Number of overlapping snps:", len(common_snps))

    print(f"New pvalue: {new_pvalue}")

    print("n snps prev - should equal number in the study: ", n_snps_prev)
    
    print("n snps after filtering: ", n_snps_filtered)
    print("sd", sd)


    print('nans in threshold')
    print(threshold_df[threshold_df['p'].isna()])

    print("nans in coding")
    print(coding_region_df[coding_region_df['p'].isna()])

    return

if __name__ == "__main__":
    if len(sys.argv) != 3:
            print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        sd = sys.argv[2]
        main(directory_gwas_combined_files, sd)

