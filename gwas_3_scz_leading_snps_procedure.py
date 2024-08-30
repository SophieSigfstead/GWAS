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

def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd):

    sd = float(sd) 

    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]
    
    threshold_set = set()
   
    non_threshold_set = set()

    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]), low_memory = False)
    non_match_list = list(first_file_df[(first_file_df['snp'].isna()) & (first_file_df['P'].notna())]['snp'])
    print("number of snps", len(first_file_df))
    print("number of snps in gwas not in 1KG: ", len(non_match_list))
    non_gwas_list = list(first_file_df[(first_file_df['P'].isna()) & (first_file_df['snp'].notna())]['snp'])
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
        chrom = row['chr']
        pos = row['pos']
        if chrom in trees:
            return bool(trees[chrom][pos])
        return False

    # Apply the function to check for overlaps
    first_file_df['in_coding_region'] = first_file_df.apply(is_within_exon, axis=1)
   
    coding_region_set = set(first_file_df[first_file_df['in_coding_region']]['snp'])
    first_file_df[first_file_df['in_coding_region']].to_csv("gwas_3_scz_intermediate_files/coding_snps.csv") '''
    
    coding_snps_df = pd.read_csv("gwas_3_scz_intermediate_files/coding_snps.csv")
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

    threshold_df = first_file_df[(first_file_df['snp'].isin(threshold_set)) & (first_file_df['snp'].notna()) & (first_file_df['P'].notna()) & (~first_file_df['snp'].isin(coding_region_set))]

    non_threshold_df = first_file_df[first_file_df['snp'].isin(non_threshold_set) &  (first_file_df['snp'].notna()) & (first_file_df['P'].notna()) & ~first_file_df['snp'].isin(coding_region_set)]
    
    # take out snps that met another threshold from non-threshold-df
    common_marker_names = set(threshold_df['snp']).intersection(set(non_threshold_df['snp']))
    if common_marker_names:
        non_threshold_df = non_threshold_df[~non_threshold_df['snp'].isin(common_marker_names)]

    non_threshold_df.to_csv(f'gwas_3_scz_intermediate_files/gwas_3_non_threshold_snps_sd={sd}.csv') 
    coding_region_df = first_file_df[(first_file_df['snp'].isin(coding_region_set)) & (first_file_df['P'].notna())]
    # Modification: the n_snps_prev should be equal to the length of the threshold and non-threshold df together
    n_snps_prev = len(threshold_df) + len(non_threshold_df) + len(coding_region_df)
    print("n snps prev - should equal number in the study: ", n_snps_prev)
    
    n_snps_filtered = len(threshold_df) + len(coding_region_df) 
    print("n snps after filtering: ", n_snps_filtered)
    
    new_pvalue = 5e-8 * (n_snps_prev / n_snps_filtered)

    print("New p-value threshold: ", new_pvalue)
    
    threshold_df['P'] = pd.to_numeric(threshold_df['P'], errors='coerce')
    print(threshold_df[threshold_df['P'].isna()])

    coding_region_df['P'] = pd.to_numeric(coding_region_df['P'],  errors='coerce')
    print(coding_region_df[coding_region_df['P'].isna()])
    #coding_region_df = coding_region_df[coding_region_df['P'] < new_pvalue]
    
    threshold_and_coding = pd.concat([threshold_df, coding_region_df], ignore_index=True)

    threshold_and_coding_sorted = threshold_and_coding.sort_values('P').reset_index()
    threshold_and_coding_sorted.to_csv(f'gwas_3_scz_intermediate_files/gwas_3_threshold_and_coding_snps_sd={sd}.csv') 
    #threshold_and_coding_sorted = pd.read_csv('gwas_3_scz_intermediate_files/gwas_3_threshold_and_coding_snps_sd=1.0.csv')

    #first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]), low_memory = False)
    #threshold_and_coding_sorted = first_file_df[(first_file_df['P'].notna()) & (first_file_df['snp'].notna())]
    plink_path = "../tools/plink2/plink2"
    test_file_path = './1KG_genotypes/all_phase3_bin_eur_no_duplicates'
   
    snp_list_file_2 = f"gwas_3_scz_intermediate_files/gwas_3_threshold_and_coding_snps_2.txt"
    snp_list_for_ld = f"gwas_3_scz_intermediate_files/gwas_3_threshold_and_coding_snps_ids.txt"
    #df_all = pd.concat(threshold_and_coding_sorted, non_threshold_df, ignore_index=True)
    df_renamed = threshold_and_coding_sorted.rename(columns={'snp': 'ID', 'P': 'P'})

    # select only the ID and P columns for the plink clumping
    df_selected = df_renamed[['ID', 'P']].reset_index()
    df_selected["P"] = df_selected["P"].astype(float)
    df_selected = df_selected[(df_selected['ID'].notna()) & (df_selected['P'].notna())]
    # select only id for ld operation
    df_selected_snps_only = df_selected[['ID']]
   
    # Save to a text file
    df_selected.to_csv(snp_list_file_2, sep='\t', index=False)
    df_selected_snps_only.to_csv(snp_list_for_ld, sep = "\t", index=False)
    
    plink_path = "../tools/plink2/plink2"
    
    # Run PLINK for clumping
    subprocess.run([
    f'{plink_path}', 
    '--bfile', f'{test_file_path}', 
    '--clump', snp_list_file_2, 
    '--clump-p1', '1e-4',
    '--clump-p2', '1e-4',
    '--clump-r2', '0.1', 
    '--clump-kb', '3000', 
    '--out', './gwas_3_scz_intermediate_files/clumped_results_2'
    ])

    # Load clumped results
    clumped_results = pd.read_csv("./gwas_3_scz_intermediate_files/clumped_results_2.clumps", delim_whitespace=True)

    filtered_df= clumped_results[((clumped_results['#CHROM']==6) & ((clumped_results['POS']>=26000000) & 
                              (clumped_results['POS'] <= 34000000)))]
    print("Filtered DF")
    print(filtered_df.head())

    
    head_snp = filtered_df.iloc[0]['ID']
    print("Head snp", head_snp)
    clumped_results_intermediate = clumped_results[((clumped_results['#CHROM']!=6)|((clumped_results['POS'] < 26000000) |
                              (clumped_results['POS'] > 34000000)))|(clumped_results['ID'] == head_snp)].sort_values(by="#CHROM")
    
    print("CLUMPED RESULTS")
    print(clumped_results_intermediate.head())
    print("LENGTH CLUMPED:", len(clumped_results_intermediate))

    clumped_results_intermediate.to_csv('./gwas_3_scz_intermediate_files/clumped_results_intermediate.csv')
    clumped_results_intermediate = clumped_results_intermediate[clumped_results_intermediate['P']<new_pvalue]

     # define the loci according to their procedure (r2 = 0.6, 500kb)

    clumped_results_with_borders_filtered = pd.DataFrame(columns = ["chr", "pos", "snp", "P", "left_border", "right_border"] )

    def define_locus(snp, chr, pos):
        ld_data = pd.read_csv(f'./gwas_3_scz_intermediate_files/ld_output_files/ld_loci_output_chr={chr}_pos={pos}.vcor', sep = '\t')
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

    for snp in list(clumped_results_intermediate["ID"]):
        row = clumped_results_intermediate[clumped_results_intermediate['ID']==snp].iloc[0]
        pos = int(row["POS"])
        chr = row['#CHROM']
        P_val = row["P"]
        start = pos - 3000000
        end = pos + 3000000

        # make the bed for each chromosome
        plink_make_bed_command = [
        f'{plink_path}', 
        '--bfile', test_file_path,  
        '--chr', str(chr),  
        '--from-bp', str(start),  
        '--to-bp', str(end),  
        '--make-bed',  
        '--out', f'./gwas_3_scz_intermediate_files/chr_bed_files/chr_bed_file_chr={chr}_pos={pos}'
        ]
        subprocess.run(plink_make_bed_command, check=True)

        # compute ld scores for each snp 
        plink_ld_command = [
            f'{plink_path}', 
            "--bfile", f'./gwas_3_scz_intermediate_files/chr_bed_files/chr_bed_file_chr={chr}_pos={pos}',  
            "--ld-window", "99999",
            "--ld-window-kb", "3000",
            "--ld-window-r2", "0.6",
            "--r2-unphased",
            "--ld-snp-list",  f'{snp_list_for_ld}',
            "--out", f"./gwas_3_scz_intermediate_files/ld_output_files/ld_loci_output_chr={chr}_pos={pos}"
        ]
        subprocess.run(plink_ld_command)

        # define locus based on ld information
        left_border, right_border = define_locus( snp, chr, pos)
        clumped_results_with_borders_filtered=clumped_results_with_borders_filtered._append({"chr":chr, "pos": pos, "snp": snp, "P": P_val, "left_border": left_border, "right_border": right_border}, ignore_index=True)

    print("Clumped Results with Borders")
    print(clumped_results_with_borders_filtered.head())
    print("Num rows:")
    print(len(clumped_results_with_borders_filtered))
    clumped_results_with_borders_filtered.to_csv(f'./gwas_3_scz_intermediate_files/clumped_results_intermediate_2_sd={sd}.csv')
    
    # perform merging of loci w.in 250kb

    def merge_loci(df):
    # Sort the DataFrame by chromosome and left_border
        df = df.sort_values(by=['chr', 'left_border']).reset_index(drop=True)
        merged_loci = []
        current_locus = df.iloc[0].copy()  # Make a copy to avoid modifying the original DataFrame
        
        for i in range(1, len(df)):
            next_locus = df.iloc[i]
            
            if (current_locus['chr'] == next_locus['chr']) and \
            (next_locus['left_border'] <= current_locus['right_border'] + 250000):
                
                print("MERGING")
                print(current_locus)
                print(next_locus)
                # Update borders
                current_locus['left_border'] = min(current_locus['left_border'], next_locus['left_border'])
                current_locus['right_border'] = max(current_locus['right_border'], next_locus['right_border'])

                
                
                # Update SNP and P value if next_locus has a lower P value
                if next_locus['P'] < current_locus['P']:
                    current_locus['snp'] = next_locus['snp']
                    current_locus['P'] = next_locus['P']
            
            else:
                # Append the merged locus to the list
                merged_loci.append(current_locus)
                current_locus = next_locus.copy()  # Start a new current locus
        
        # Append the last locus
        merged_loci.append(current_locus)
        
        # Convert the list of merged loci back to a DataFrame
        merged_df = pd.DataFrame(merged_loci)
        
        return merged_df

    clumped_results_with_merged_borders = merge_loci(clumped_results_with_borders_filtered)

    print((clumped_results_with_merged_borders).head(10))

    clumped_results_with_merged_borders.to_csv(f'./gwas_3_scz_result_files/filtered_snps_sd={sd}.csv')
    
    ground_truth_snps = pd.read_csv('./gwas_3_scz_result_files/filtered_snps_sd=0.0.csv')

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
                    
        return overlap_count

    print("number of snps found with our method", len(clumped_results_with_merged_borders))

    filtered_in_coding_region = clumped_results_with_merged_borders[clumped_results_with_merged_borders['snp'].isin(coding_region_set)]

    print("number of coding snps" , len(filtered_in_coding_region))
    print(filtered_in_coding_region)
    filtered_in_coding_region.to_csv(f'gwas_3_scz_intermediate_files/our_coding_snps_sd={sd}')

    ground_truth_coding = ground_truth_snps[ground_truth_snps['snp'].isin(coding_region_set)]
    print("number of coding snps originally in ground truth: ", len(ground_truth_coding))
    print(ground_truth_coding)
    ground_truth_coding.to_csv(f'gwas_3_scz_intermediate_files/ground_truth_coding_snps_sd={sd}')

    overlap_count = count_overlaps(ground_truth_snps, clumped_results_with_merged_borders)

    print(f"Number of overlapping loci: {overlap_count}")

    common_snps = set(ground_truth_snps['snp']).intersection(clumped_results_with_merged_borders['snp'])

    print(f"Number of overlapping snps:", len(common_snps))

    print(f"New pvalue: {new_pvalue}")

    print("n snps prev - should equal number in the study: ", n_snps_prev)
    
    print("n snps after filtering: ", n_snps_filtered)
    print("sd", sd)


    print('nans in threshold')
    print(threshold_df[threshold_df['P'].isna()])

    print("nans in coding")
    print(coding_region_df[coding_region_df['P'].isna()])

    return

if __name__ == "__main__":
    if len(sys.argv) != 4:
            print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        directory_gwas_sig_snps = sys.argv[2]
        sd = sys.argv[3]
        main(directory_gwas_combined_files, directory_gwas_sig_snps, sd)
      
   
