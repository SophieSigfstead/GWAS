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
    # Ensure the dataframe is sorted by the left border
    df = df.sort_values(by='left_border').reset_index(drop=True)
    
    # Initialize a list to hold the merged intervals
    merged_intervals = []
    
    # Start with the first interval
    current_interval = df.iloc[0][['chr', 'left_border', 'right_border']].tolist()
    
    # Iterate over the remaining intervals
    for i in range(1, len(df)):
        chr_value, left_border, right_border = df.iloc[i][['chr', 'left_border', 'right_border']]
        
        # Check if the current interval overlaps with the next one
        if left_border <= current_interval[2]:
            # If they overlap, merge the intervals by updating the right border
            current_interval[2] = max(current_interval[2], right_border)
        else:
            # If they do not overlap, add the current interval to the list
            merged_intervals.append(current_interval)
            # Start a new current interval
            current_interval = [chr_value, left_border, right_border]
    
    # Add the last interval
    merged_intervals.append(current_interval)
    
    # Convert the list of merged intervals back to a DataFrame
    merged_df = pd.DataFrame(merged_intervals, columns=['chr', 'left_border', 'right_border'])
    
    return merged_df


def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, gcta_path):
    sd = float(sd) 
    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))
    our_snps_list = list(first_file_df['snp'])
    their_snps = pd.read_csv('./gwas_2_intermediate_files/their_loci_list.csv')
    their_snps_list  = list(their_snps['snp'])
    snps_not_in_our_list = result = [item for item in their_snps_list if item not in our_snps_list]
    print("number of snps missing from our list:")
    print(len(snps_not_in_our_list))
   

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
    if sd == 0:
        new_pvalue = 5e-8
    print("New p-value threshold: ", new_pvalue)
    print("len(threshold_df) before filtering:",  len(threshold_df))

    

    threshold_df = threshold_df[(threshold_df['P'].notna()) & (threshold_df['P'] < new_pvalue)]
    print("len(threshold_df) after filtering:",  len(threshold_df))
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
    threshold_and_coding = threshold_and_coding.rename(columns={'p': 'P'})
    threshold_and_coding = threshold_and_coding.rename(columns={'SNP': 'ID'})
    threshold_and_coding_extract= threshold_and_coding[['ID', 'P']]
    threshold_and_coding_extract.to_csv('./gwas_2_intermediate_files/significant_variants.csv', sep='\t', index=False)
    threshold_and_coding = threshold_and_coding.rename(columns={'ID': 'SNP'})
    test_file_path = './1KG_genotypes/all_phase3_bin'
    result_df = pd.DataFrame()

    if sd == 0:
        new_pvalue=5e-8 


    plink_path = "../tools/plink2/plink2"
    # Run PLINK for clumping
    subprocess.run([
    f'{plink_path}', 
    '--bfile', f'{test_file_path}', 
    '--clump', 'gwas_2_intermediate_files/significant_variants.csv', 
    '--clump-p1', '1e-4',
    '--clump-p2', '1e-4',
    '--clump-r2', '0.1', 
    '--clump-kb', '3000', 
    '--out', './gwas_2_intermediate_files/clumped_results'
    ])

    # Load clumped results
    clumped_results = pd.read_csv("./gwas_2_intermediate_files/clumped_results.clumps", delim_whitespace=True)
    print(clumped_results.head())

    test_file_path = './1KG_genotypes/all_phase3_bin'

    os.makedirs('./gwas_2_intermediate_files/chromosome_data_files', exist_ok=True)

    os.makedirs("./gwas_2_intermediate_files/gcta_output_directory", exist_ok=True)

    loci_data = pd.DataFrame(columns = ['chr', 'left_border', 'right_border'])

    for index, row in clumped_results.iterrows():
        locus_snps = row['SP2'].split(',')  # Extract SNPs in the locus
        locus_data = threshold_and_coding[threshold_and_coding['SNP'].isin(locus_snps)]
        if not locus_data.empty:
        #print(locus_data.head())
            left_border = int(locus_data['pos'].min())
            right_border = int(locus_data['pos'].max())
        else:
            left_border = int(row['POS'])
            right_border = int(row['POS'])
        chr = row['#CHROM']
        loci_data = loci_data._append({ 'chr': chr, 'left_border': left_border, 'right_border': right_border}, ignore_index=True)

    loci_data = loci_data.sort_values(by=['chr', 'left_border'])
    print(loci_data)

    # merge the loci in loci_data
    print("len loci data", len(loci_data))
    merged_loci = pd.DataFrame(columns = ['chr', 'left_border', 'right_border'])
    
    for chr in range(1, 23):  
        chr_loci = loci_data[loci_data['chr']==chr]
        if not chr_loci.empty: # provides dataframe with chr (same for all rows), left border, right border
            chr_loci_merged = merge_intervals(chr_loci) # must return dataframe with chr, left border, right border
            merged_loci = pd.concat([merged_loci, chr_loci_merged], ignore_index=True)

    print("len merged loci", len(merged_loci))
    merged_loci.to_csv('./merged_loci.csv')
    print(merged_loci)

    for _,row in merged_loci.iterrows():
        chr = row['chr']
        left_border = row['left_border']
        right_border = row['right_border']
        filtered_file_path = f'./gwas_2_intermediate_files/chromosome_data_files/snp_list_chr={chr}_left_border={left_border}.ma'
        filtered_file_df = threshold_and_coding[(threshold_and_coding['chr'] == chr) & (threshold_and_coding['pos'] >= left_border) & (threshold_and_coding['pos'] <= right_border) ]
        filtered_file_df.to_csv(filtered_file_path,sep='\t', index=False)
        gcta_output_file = os.path.join("./gwas_2_intermediate_files/gcta_output_directory", f'test_chr={chr}_left_border={left_border}')

        gcta_command = f'{gcta_path} --bfile {test_file_path} --chr {chr} --cojo-file {filtered_file_path} --cojo-slct --out {gcta_output_file}'
        subprocess.run(gcta_command, shell=True)

        gcta_result_file = f'{gcta_output_file}.jma.cojo'
        if os.path.exists(gcta_result_file):
            chr_result_df = pd.read_csv(gcta_result_file, sep='\t')
            result_df = pd.concat([result_df, chr_result_df], ignore_index=True)
        
    combined_results_path = os.path.join("./gwas_2_intermediate_files/gcta_output_directory", 'combined_results.ma')
    #result_df.to_csv(combined_results_path, sep='\t', index=False)

    result_df = pd.read_csv(combined_results_path,  sep ='\t')
    merged_loci = pd.read_csv('./merged_loci.csv')
    print(merged_loci.head())

    print('number of snps found', len(result_df))
    print(result_df.head())
    common_snps = set(result_df['SNP']).intersection(set(their_snps_list))
    print('number of snps in common', len(list(common_snps)))
    
    # Function to check overlaps
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

    their_loci_list = pd.read_csv('./gwas_2_intermediate_files/their_loci_list.csv')
    print(their_loci_list.head())
    # Count overlaps between the two DataFrames
    overlap_count = count_overlaps(merged_loci, their_loci_list)


    print(f"Number of overlapping intervals: {overlap_count}")
    #print('Number of common SNPs:', len(common_snps))
    #print("New p-value threshold: ", new_pvalue)
    #print(common_snps)
    
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
