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


def define_ld_boundaries(ld_df, snp):
    subset = ld_df[(ld_df['SNP_A'] == snp) | (ld_df['SNP_B'] == snp)]
    if subset.empty:
        return None, None
    left_boundary = min(subset['BP_A'].min(), subset['BP_B'].min())
    right_boundary = max(subset['BP_A'].max(), subset['BP_B'].max())
    return left_boundary, right_boundary

def merge_intervals(intervals):
    tree = IntervalTree(Interval(begin, end) for begin, end in intervals)
    merged = sorted(tree.merge_overlaps())
    return [(interval.begin, interval.end) for interval in merged]


def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, gcta_path):
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

    threshold_df = threshold_df[(threshold_df['P'].notna()) & (threshold_df['P'] < new_pvalue)]
    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region']) & (first_file_df['P'] < new_pvalue)]
    
    threshold_and_coding = pd.concat([threshold_df, coding_region_df], ignore_index=True)

    # reformat the threshold and coding SNPs list into the required file format for cojo
    # columns must be SNP A1 A2 freq b se p N  and format is .ma
    print("threshold and coding head:")
    print(threshold_and_coding.head())

    threshold_and_coding = threshold_and_coding.drop(columns=['A1', 'A2', 'SAD411', 'in_coding_region', ])

    column_mapping = {
    'ref': 'A1',
    'alt': 'A2',
    'snp': 'SNP',
    'P': 'p',
    'StdErrLogOR': 'se',
    'LogOR': 'b',
    'Freq': 'freq',
    'MarkerName':'MarkerName'
    }

    threshold_and_coding.rename(columns=column_mapping, inplace=True)
    threshold_and_coding['N'] = 807553
    desired_order = ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N', 'chr', 'pos', "MarkerName"]

    threshold_and_coding = threshold_and_coding[desired_order]

     # sample size

    result_df = pd.DataFrame()
    output_directory_for_gcta_files = './gcta_output_directory'
    os.makedirs(output_directory_for_gcta_files, exist_ok=True) 

    for chr_num in range(1, 23):
        filtered_snp_list = threshold_and_coding[threshold_and_coding['chr'] == chr_num]
        filtered_file_path = os.path.join(output_directory_for_gcta_files, f'snp_list_chr{chr_num}.ma')
        filtered_snp_list.to_csv(filtered_file_path, sep='\t', index=False)
        
        gcta_output_file = os.path.join(output_directory_for_gcta_files, f'test_chr{chr_num}')

        test_file_path = './gcta-1.94.1-linux-kernel-3-x86_64/test.fam'

        gcta_command = f'{gcta_path} --bfile {test_file_path} --chr {chr_num} --maf 0.01 --cojo-file {filtered_file_path} --cojo-slct --out {gcta_output_file}'
        subprocess.run(gcta_command, shell=True)
        
        gcta_result_file = f'{gcta_output_file}.jma.cojo'
        if os.path.exists(gcta_result_file):
            chr_result_df = pd.read_csv(gcta_result_file, sep='\t')
            result_df = pd.concat([result_df, chr_result_df], ignore_index=True)
        
        print(f'Processed chromosome {chr_num} and saved results to {gcta_output_file}')

    result_df.to_csv(os.path.join(output_directory_for_gcta_files, 'combined_results.ma'), sep='\t', index=False)
    print("Processing complete for all chromosomes.")
    
    os.makedirs('./filtered_snps_gwas_2', exist_ok=True) 
    result_df.to_csv(f'./filtered_snps_gwas_2/filtered_snps_gwas_2_threshold={sd}.csv', index=False)

    print('number of snps found', len(result_df))
    
    # GWAS 2
    # compare SNPs with the significant SNPs file
    sig_snps_df = pd.read_csv(directory_gwas_sig_snps)
    
    sig_snps_set = set(sig_snps_df['MarkerName'])
    
    result_snps_set = set(result_df['SNP'])


    common_snps = result_snps_set.intersection(sig_snps_set)
    print('Number of common SNPs:', len(common_snps))

    #command = f"./plink-1.07-x86_64/plink --bfile ./<DATA> --r2 --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.8 --out ./ld_results"
    # plink --bfile your_data --r2 --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.8 --out ld_results
    #subprocess.run(command, shell=True)
    
    ld_file = './ld_results.ld'  # You should calculate this file with PLINK as described earlier
    ld_results = pd.read_csv(ld_file, delim_whitespace=True)
    
    intervals = []
    for snp in result_snps_set:
        left, right = define_ld_boundaries(ld_results, snp)
        if left is not None and right is not None:
            intervals.append((left, right))
    
    merged_intervals = merge_intervals(intervals)
    
    other_snps_positions = sig_snps_df.set_index('snp')['pos'].to_dict()
    
    common_loci = []
    for snp, pos in other_snps_positions.items():
        for left, right in merged_intervals:
            if left <= pos <= right:
                common_loci.append(snp)
                break
    
    print('Number of common loci:', len(set(common_loci)))
    
    return
    
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd> < directory_for_gcta64>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        directory_gwas_sig_snps = sys.argv[2]
        sd = sys.argv[3]
        directory_for_gcta64 = [4]
        main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, directory_for_gcta64)
    sys.exit()
