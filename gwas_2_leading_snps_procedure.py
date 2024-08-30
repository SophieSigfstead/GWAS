import pandas as pd
import os
import sys
from intervaltree import IntervalTree, Interval
import subprocess

def threshold_calculator(df, sd, track):
    df_with_sad = df[(df[track].notna()) & (df[track].abs() < 1) & (df[track].abs() >= 0)]
  
    std_dev = df_with_sad[track].std()
    mean = df_with_sad[track].mean()

    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)

    return upper_threshold, lower_threshold


def merge_intervals(df):
    merged = pd.DataFrame(columns=['Chr', 'left_border', 'right_border'])
    
    # Group by chromosome
    for chr, group in df.groupby('Chr'):
        group = group.sort_values('left_border')
        start = None
        end = None
        for _, row in group.iterrows():
            if start is None:
                start = row['left_border']
                end = row['right_border']
            elif row['left_border'] <= end:  # Overlapping or adjacent intervals
                end = max(end, row['right_border'])
            else:
                merged = merged._append({'Chr': chr, 'left_border': start, 'right_border': end}, ignore_index=True)
                start = row['left_border']
                end = row['right_border']
        # Append the last interval
        if start is not None:
            merged = merged._append({'Chr': chr, 'left_border': start, 'right_border': end}, ignore_index=True)
    
    return merged

def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, gcta_path):
    sd = float(sd)
    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]
    
    threshold_set = set()
    non_threshold_set = set()

    exon_regions = pd.read_csv('./genome_assembly/exon_regions_v2.csv')
    
    # Read the first file to initialize the in_coding_region column
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))

    # create interval trees for each exon in exon regions 
    trees = {}
    for _, row in exon_regions.iterrows():
        chr = row['chr']
        if chr not in trees:
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

    
    # Modification: the n_snps_prev should be equal to the length of the threshold and non-threshold df together
    n_snps_prev = len(threshold_df) + len(non_threshold_df) + len(coding_region_set)
    print("n snps prev", n_snps_prev)

    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region'])]
    
    n_snps_filtered = len(threshold_df) + len(coding_region_df)
    print("n snps filtered", n_snps_filtered)
    
    new_pvalue = 5e-8 * (n_snps_prev / n_snps_filtered)
    print("New p-value threshold: ", new_pvalue)

    threshold_df = threshold_df[(threshold_df['P'].notna()) & (threshold_df['P'] < new_pvalue)]
    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region']) & (first_file_df['P'] < new_pvalue)]
    threshold_df_old = threshold_df.copy()
    coding_region_df_old = coding_region_df.copy()
    coding_region_df_old.to_csv('./gwas_2_intermediate_files/coding_regions_df.csv')
    threshold_and_coding_old = pd.concat([threshold_df_old, coding_region_df_old], ignore_index = True)
    
    threshold_and_coding = pd.concat([threshold_df, coding_region_df], ignore_index=True)

    # reformat the threshold and coding SNPs list into the required file format for cojo
    # columns must be SNP A1 A2 freq b se p N  and format is .ma
    print("threshold and coding head:")
    print(threshold_and_coding.head())

    threshold_and_coding = threshold_and_coding.drop(columns=['A1', 'A2', 'SAD411'])

    column_mapping = {
    'ref': 'A1',
    'alt': 'A2',
    'snp': 'SNP',
    'P': 'p',
    'StdErrLogOR': 'se',
    'LogOR': 'b',
    'Freq': 'freq',
    'MarkerName':'MarkerName',
    'in_coding_region': 'in_coding_region'
    }

    threshold_and_coding.rename(columns=column_mapping, inplace=True)
    threshold_and_coding['N'] = 807553
    desired_order = ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N', 'chr', 'pos', "MarkerName", 'in_coding_region']

    threshold_and_coding = threshold_and_coding[desired_order]
    

    threshold_and_coding.to_csv('threshold_and_coding.csv') 
    non_threshold_df.to_csv('non_threshold.csv') 
    threshold_and_coding_old.to_csv('threshold_and_coding_old.csv') 
    
    threshold_and_coding_old = pd.read_csv('threshold_and_coding_old.csv')
    threshold_and_coding = pd.read_csv('threshold_and_coding.csv')
    non_threshold_df= pd.read_csv('non_threshold.csv')
    desired_order = ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N', 'chr', 'pos', "MarkerName",'in_coding_region' ]
    threshold_and_coding = threshold_and_coding[desired_order]
    print(threshold_and_coding.head())

    test_file_path = './1KG_genotypes/all_phase3_bin'
    result_df = pd.DataFrame()
    output_directory_for_gcta_files = './gcta_output_directory'
    os.makedirs(output_directory_for_gcta_files, exist_ok=True) 

    for chr_num in range(1, 23):
        filtered_snp_list = threshold_and_coding[threshold_and_coding['chr'] == chr_num]
        if len(filtered_snp_list) > 0:
            filtered_file_path = os.path.join(output_directory_for_gcta_files, f'snp_list_chr{chr_num}.ma')
            filtered_snp_list.to_csv(filtered_file_path, sep='\t', index=False)
            
            gcta_output_file = os.path.join(output_directory_for_gcta_files, f'test_chr{chr_num}')

            gcta_command = f'{gcta_path} --bfile {test_file_path} --chr {chr_num} --cojo-file {filtered_file_path} --cojo-slct --out {gcta_output_file}'
            subprocess.run(gcta_command, shell=True)
            
            gcta_result_file = f'{gcta_output_file}.jma.cojo'
            if os.path.exists(gcta_result_file):
                chr_result_df = pd.read_csv(gcta_result_file, sep='\t')
                result_df = pd.concat([result_df, chr_result_df], ignore_index=True)
            
            #print(f'Processed chromosome {chr_num} and saved results to {gcta_output_file}')
    combined_results_path = os.path.join(output_directory_for_gcta_files, 'combined_results.ma')
    result_df.to_csv(combined_results_path, sep='\t', index=False)
    print("Processing complete for all chromosomes.")
    
    os.makedirs('./filtered_snps_gwas_2', exist_ok=True) 
    result_df.to_csv(f'./filtered_snps_gwas_2/filtered_snps_gwas_2_threshold={sd}.csv', index=False)

    print('number of snps found', len(result_df))
    print(result_df.head())
    
    
    # GWAS 2
    # compare SNPs with the significant SNPs file
    sig_snps_df = pd.read_csv(directory_gwas_sig_snps)
    
    sig_snps_set = set(sig_snps_df['MarkerName'])
    sig_snps_list = list(sig_snps_set)
  
    result_snps_set = set(result_df['SNP'])

    result_snps_list = list(result_snps_set)

    sig_snps_list = list(sig_snps_set)

    # Save the list to a text file, each SNP on a new line
    with open('./gwas_2_intermediate_files/gcta_result_snps_list.txt', 'w') as file:
        for snp in result_snps_list:
            file.write(f"{snp}\n") 


    all_df = pd.concat([threshold_and_coding_old, non_threshold_df], ignore_index = True)

    file_path = './gwas_2_intermediate_files/gcta_result_snps_list.txt'

    
    result_snp_list = []

    
    with open(file_path, 'r') as file:
        for line in file:
          
            snp = line.strip()
            
            if snp:
                result_snp_list.append(snp)
    
    all_df.to_csv("./gcta_output_directory/all_snps_df.csv")
    common_snps = set(result_snp_list).intersection(sig_snps_set)
    print('Number of common SNPs:', len(common_snps))

    plink_path = "../tools/plink2/plink2"

    

    LD_THRESHOLD = 0.75

    our_loci_list = pd.DataFrame(columns = ['snp', 'chr', 'pos', 'left_border', 'right_border', 'in_coding_region'])

    for snp in result_snp_list:
        row = all_df[all_df['snp'] == snp].iloc[0]
        chr = int(row['chr'])
        pos = row['pos']
        in_coding_region = row['in_coding_region']
        start = pos - 1500000
        end = pos + 1500000
        ld_matrix_command = f'{plink_path} --bfile {test_file_path} --chr {chr} --from-bp {start} --to-bp {end} --make-bed --out ./gcta_output_directory/region_data'
        subprocess.run(ld_matrix_command, shell=True)
        region_path = './gcta_output_directory/region_data'
        r2_matrix_command = f'{plink_path} --bfile {region_path} --r2-unphased --out ./gcta_output_directory/ld_matrix'
        subprocess.run(r2_matrix_command, shell=True)
        ld_matrix_path = './gcta_output_directory/ld_matrix.vcor'
        ld_matrix = pd.read_csv(ld_matrix_path, delim_whitespace=True)
        high_ld_with_target = ld_matrix[(ld_matrix['ID_A'] == snp) & (ld_matrix['UNPHASED_R2'] > LD_THRESHOLD)]
        print(high_ld_with_target.head())
        if (len(high_ld_with_target)>0):
            pos_right = high_ld_with_target['POS_B'].astype(int).max()
            pos_left = high_ld_with_target['POS_B'].astype(int).min()
        else:
            pos_right = pos
            pos_left = pos
        our_loci_list = our_loci_list._append({'snp': snp, 'chr': chr, 'pos': pos, 'left_border': pos_left, 'right_border': pos_right,'in_coding_region': in_coding_region}, ignore_index=True)

    print(our_loci_list.head())
    print("number of coding snps:")
    print(len(our_loci_list[our_loci_list['in_coding_region']==True]))
    our_loci_list.to_csv(f'filtered_snps_gwas_2/gwas_2_loci_lists/out_loci_list_sd={sd}.csv')


    their_loci_list = pd.DataFrame(columns = ['snp', 'chr', 'pos', 'left_border', 'right_border', 'in_coding_region'])
    snp_info_df = pd.read_csv('./gwas_2_intermediate_files/snp_info_df')

    for snp in sig_snps_list:
        row = snp_info_df[snp_info_df['SNP'] == snp].iloc[0]
        chr = int(row['Chr'])
        pos = row['BP']
        if len(all_df[all_df['snp'] == snp] !=0):
            in_coding_region = all_df[all_df['snp'] == snp]['in_coding_region'].iloc[0].astype(bool)
        else:
            print("Could not find")
            in_coding_region = False
        start = pos - 1500000
        end = pos + 1500000
        ld_matrix_command = f'{plink_path} --bfile {test_file_path} --chr {chr} --from-bp {start} --to-bp {end} --make-bed --out ./gcta_output_directory/region_data'
        subprocess.run(ld_matrix_command, shell=True)
        region_path = './gcta_output_directory/region_data'
        r2_matrix_command = f'{plink_path} --bfile {region_path} --r2-unphased --out ./gcta_output_directory/ld_matrix'
        subprocess.run(r2_matrix_command, shell=True)
        ld_matrix_path = './gcta_output_directory/ld_matrix.vcor'
        ld_matrix = pd.read_csv(ld_matrix_path, delim_whitespace=True)
        high_ld_with_target = ld_matrix[(ld_matrix['ID_A'] == snp) & (ld_matrix['UNPHASED_R2'] > LD_THRESHOLD)]
        if (len(high_ld_with_target)>0):
            pos_right = high_ld_with_target['POS_B'].astype(int).max()
            pos_left = high_ld_with_target['POS_B'].astype(int).min()
        else:
            pos_right = pos
            pos_left = pos
        their_loci_list = their_loci_list._append({'snp': snp, 'chr': chr, 'pos': pos, 'left_border': pos_left, 'right_border': pos_right,'in_coding_region': in_coding_region}, ignore_index=True)

    print(their_loci_list.head())
    print("number of coding snps in their loci list:")
    print(len(their_loci_list[their_loci_list['in_coding_region'] == True]))
  
    their_loci_list.to_csv(f'filtered_snps_gwas_2/gwas_2_loci_lists/their_loci_list.csv')
    
    # Function to check overlaps
    def count_overlaps(df1, df2):
        
        overlaps = 0
        
        for i, row1 in df1.iterrows():
            left1, right1 = row1['left_border'], row1['right_border']
            chr1  = row1['chr']
        
            # Check for overlap with each interval in df2
            for j, row2 in df2.iterrows():
                left2, right2 = row2['left_border'], row2['right_border']
                chr2  = row2['chr']
            
                # Check if the intervals overlap
                if (left1 <= right2 and right1 >= left2) and (chr1 == chr2):
                    overlaps += 1
                    break
        
        return overlaps

    # Count overlaps between the two DataFrames
    overlap_count = count_overlaps(our_loci_list, their_loci_list)


    print(f"Number of overlapping intervals: {overlap_count}")
    print('Number of common SNPs:', len(common_snps))
    print("New p-value threshold: ", new_pvalue)
    print(common_snps)
    
    return

    
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd> <directory_for_gcta64>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        directory_gwas_sig_snps = sys.argv[2]
        sd = sys.argv[3]
        gcta_path = sys.argv[4]
        main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, gcta_path)
    sys.exit()
