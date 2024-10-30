# Purpose of file: To perform an analysis  to compare the performance of GWAS 1 single tracks to random 33-track sets for SNP selection.and
# Author: Sophie Sigfstead
# 

# Purpose of file: To perform an analysis to compare the performance of GWAS 1 single tracks to random 33-track sets for SNP selection.
# Author: Sophie Sigfstead

import random
import subprocess
import os
import pandas as pd 
import matplotlib.pyplot as plt
import shutil

def compute_non_overlapping_loci(single_track_df, reference_loci):
    overlap_loci = pd.DataFrame(columns=['snp', 'chr', 'left_border', 'right_border'])
    single_track_df = single_track_df[['snp', 'chr', 'left_border', 'right_border']]
    for _, row1 in single_track_df.iterrows():
        left1, right1, chr1 = row1['left_border'], row1['right_border'], row1['chr']
        for _, row2 in reference_loci.iterrows():
            left2, right2, chr2 = row2['left_border'], row2['right_border'], row2['chr']
            if (left1 <= right2 and right1 >= left2) and (chr1 == chr2):
                overlap_loci = pd.concat([overlap_loci, pd.DataFrame([row2])], ignore_index=True)
                break
    missing_loci = reference_loci[~(reference_loci['snp'].isin(overlap_loci['snp']))]
    return overlap_loci, missing_loci

def clear_directory(directory_path):
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)
    os.makedirs(directory_path)

def main():
    # Step 1. Define a dictionary that houses random track sets to run the single track analysis on. 
    # constants to set up the dictionary
    num_track_sets = 2
    num_tracks_per_set = 33
    track_set_dictionary = {}

    # brain tracks - to be excluded from the random track sets
    brain_tracks = {0, 1, 9, 76, 78, 80, 81, 172, 179, 216, 240, 261, 278,
                    319, 326, 338, 355, 370, 403, 411, 421, 458, 462, 469,
                    499, 524, 552, 580, 582, 602, 644, 669}
    
    # Create the keys - these will be the ids of the track sets
    keys = [chr(i) for i in range(ord('a'), ord('a')+num_track_sets)]

    # Create the list of tracks to choose from (non-brain)
    valid_numbers = [i for i in range(684) if i not in brain_tracks]

    for key in keys:
        track_set_dictionary[key] = random.sample(valid_numbers, num_tracks_per_set)

    # to hold the results of each track set
    directories = []

    for key, value in track_set_dictionary.items():
        # Step 1. Match up the summary statistics with the data from 1000 genomes
        command = [
            'python', 
            '../gwas_1_matching.py', 
            '../gwas_1_and_2_summary_statistics_data/GCST90277450.tsv', 
            '../1000genomes_as_csv', 
            f"{value}"
        ]
        subprocess.run(command)

        # Step 2, using the matched files, run the GWAS 1 single track method
        matched_files = './gwas_random_tracks_matching'
        sd_list = [1.0]
        for sd in sd_list: 
            command = [
                'python', 
                'gwas_1_leading_SNPs_by_track.py', 
                f"{matched_files}", 
                '../gwas_1_intermediate_files/mdd_sig_snps.csv', 
                f"{sd}", 
                f"{key}"
            ]
            subprocess.run(command)
        
        directories.append(f"./gwas_1_leading_SNPs_by_track_random_results/gwas_1_leading_SNPs_by_track_random_{key}/sd={sd}")
    print(directories)
    # create the track labels info
    track_labels = "../target_dnase_ataq_tracks_labelled.csv"
    track_labels_df = pd.read_csv(track_labels)
    track_labels_df['track_col'] = 'SAD' + track_labels_df.index.astype(str)
    track_labels_df = track_labels_df[['track_col', 'target_labels', 'tissue']]

    # create list of coding snps
    coding_snps_list = "../gwas_3_scz_intermediate_files/coding_snps.csv"
    coding_snps_df = pd.read_csv(coding_snps_list)[['chr', 'pos', 'snp']]
    coding_snps_set = set(coding_snps_df['snp'])
    coding_snps_ref = {'rs9074', 'rs9262132', 'rs9586'}

    # create original snp list from the study 
    original_snps_path = '../filtered_snps_gwas_1/filtered_snps_gwas_1_threshold=0.0.csv'
    original_snps = pd.read_csv(original_snps_path)
    original_snps['left_border'] = original_snps['pos'] - 1000000
    original_snps['right_border'] = original_snps['pos'] + 1000000
    reference_loci = original_snps[['snp', 'chr', 'left_border', 'right_border']]
    original_snps_list = list(original_snps['snp'])

    # path for results from running our method with 1.0 sd on brain track set 
    gwas_1_all_track_results_path = '../filtered_snps_gwas_1/filtered_snps_gwas_1_threshold=1.0.csv'

    # results directory to save results to from the below analyses
    results_directory = f'./gwas_1_random_track_set_analysis_results/gwas_1_random_track_set_analysis_results_{key}/sd={sd}'
    clear_directory(results_directory)

    for directory in directories:
        print("Processing directory:", directory)

        # list all files in the directory 
        files = [file for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]

        # read the overlap csv files - represent the overlapping snps with the original study
        overlap_dfs = [file for file in files if file.endswith('overlap_df.csv')]
        overlap_df_list = []
        for file in overlap_dfs:
            df = pd.read_csv(os.path.join(directory, file))
            overlap_df_list.append(df)

        # Concat all overlapping datframes
        full_overlap_dfs = pd.concat(overlap_df_list, ignore_index=True)
        full_overlap_dfs = pd.merge(full_overlap_dfs, track_labels_df, on='track_col', how='inner')

        # Add summary row and save
        full_overlap_dfs = full_overlap_dfs._append(
        {'track_col': 'ALL_brain', 'matching_snps': 15, 'matching_loci': 28, 'new_coding_snps': 3, 
         'target_labels': 'all_brain', 'tissue': 'random'}, ignore_index=True).sort_values(by=['matching_loci', 'matching_snps'], ascending=False)
        full_overlap_dfs.to_csv(f'{results_directory}/full_overlap_df.csv', index=False)
    
        #read in GWAS results
        gwas_1_single_track_results = [file for file in files if file.endswith(f'sd={sd}.csv')]
        gwas_1_all_track_results = pd.read_csv("../filtered_snps_gwas_1/filtered_snps_gwas_1_threshold=1.0.csv")

        all_track_df_coding = gwas_1_all_track_results[gwas_1_all_track_results['snp'].isin(coding_snps_set)]
        non_coding_all_tracks = gwas_1_all_track_results[~gwas_1_all_track_results['snp'].isin(coding_snps_set)]

        coding_snps_comparison_df=pd.DataFrame()
        # coding comparison
        for file in gwas_1_single_track_results:
            single_track_df = pd.read_csv(os.path.join(directory, file))
            single_track_df_coding = single_track_df[single_track_df['snp'].isin(coding_snps_set)]
            track_col = [col for col in single_track_df.columns if col.startswith('SAD')][0] 
            coding_snps = set(single_track_df_coding['snp'])
            inter = len(coding_snps.intersection(coding_snps_ref))
            unique_coding_snps = coding_snps - coding_snps_ref
            coding_snps_comparison_df = coding_snps_comparison_df._append({
                'track_col': track_col, 
                'track_coding_snps': len(single_track_df_coding),
                'coding_snps': coding_snps, 
                'intersection_with_all_track_coding_snps': inter,
                'number_of_unique_snps': len(unique_coding_snps),
                'unique_coding_snps': unique_coding_snps
            }, ignore_index=True)
        

        coding_snps_comparison_df = pd.merge(coding_snps_comparison_df, track_labels_df, on='track_col', how='left')

        coding_snps_comparison_df.to_csv(f'{results_directory}/coding_snps_comparison.csv', index=False)


        # Non-coding comparison dataframe
        non_coding_comparison = pd.DataFrame(columns=[
        'track_col', 'num_snps_total', 'num_non_coding', 'num_coding', 
        'num_non_coding_intersect', 'non_coding_intersect_with_all_tracks', 
        'num_non_coding_intersect_original', 'non_coding_intersect_with_original', 
        'num_overlap_loci', 'overlap_loci', 'num_missing_loci', 'missing_loci'
        ])

        for file in gwas_1_single_track_results:
            single_track_df = pd.read_csv(os.path.join(directory, file))
            total_found = len(single_track_df)
            single_track_df_coding = single_track_df[single_track_df['snp'].isin(coding_snps_set)]
            single_track_df_non_coding = single_track_df[~single_track_df['snp'].isin(coding_snps_set)]
            
            # Track column
            track_col = [col for col in single_track_df.columns if col.startswith('SAD')][0] 
            
            # Intersections
            non_coding_intersect = set(single_track_df_non_coding['snp']).intersection(set(non_coding_all_tracks['snp']))
            non_coding_intersect_original_study = set(single_track_df_non_coding['snp']).intersection(set(original_snps_list))
            
            # Overlap and missing loci
            overlap_loci_df, missing_loci_df = compute_non_overlapping_loci(single_track_df, reference_loci)
            
            # Append data
            non_coding_comparison = non_coding_comparison._append({
                'track_col': track_col,
                'num_snps_total': total_found,
                'num_non_coding': len(single_track_df_non_coding),
                'num_coding': len(single_track_df_coding),
                'num_non_coding_intersect': len(non_coding_intersect),
                'non_coding_intersect_with_all_tracks': non_coding_intersect,
                'num_non_coding_intersect_original': len(non_coding_intersect_original_study),
                'non_coding_intersect_with_original': non_coding_intersect_original_study,
                'num_overlap_loci': len(overlap_loci_df),
                'overlap_loci': set(overlap_loci_df['snp']),
                'num_missing_loci': len(missing_loci_df),
                'missing_loci': set(missing_loci_df['snp'])
            }, ignore_index=True)

        non_coding_comparison = pd.merge(non_coding_comparison, track_labels_df, on='track_col', how='left')

        # Save non-coding comparison results
        non_coding_comparison.to_csv(f'{results_directory}/non_coding_comparison.csv', index=False)

        correlation_coding = non_coding_comparison['num_snps_total'].corr(non_coding_comparison['num_coding'])
        print(f"Correlation coefficient between total snps and number of coding snps found: {correlation_coding}\n")

        correlation_non_coding = non_coding_comparison['num_snps_total'].corr(non_coding_comparison['num_non_coding'])
        print(f"Correlation coefficient between total snps and number of non-coding snps found: {correlation_non_coding}")

        # mean_matching_loci calculation
        mean_matching_loci = full_overlap_dfs['matching_loci'].mean()
        median_matching_loci = full_overlap_dfs['matching_loci'].median()
        standard_dev_matching_loci = full_overlap_dfs['matching_loci'].std()
        print(f"Mean matching loci for {directory}: {mean_matching_loci}")
        stats_df = pd.DataFrame({
            'mean_matching_loci': [mean_matching_loci],
            'median_matching_loci': [median_matching_loci],
            'standard_dev_matching_loci': [standard_dev_matching_loci],
            'correlation_non_coding_to_total_snps_found': [correlation_non_coding],
            'correlation_coding_to_total_snps_found': [correlation_coding]
        })

        # Save to CSV in the results directory
        stats_df.to_csv(os.path.join(results_directory, 'matching_loci_stats.csv'), index=False)

        plt.figure(figsize=(8, 6))
        plt.scatter(non_coding_comparison['num_snps_total'], non_coding_comparison['num_coding'], color='b', alpha=0.6)
        plt.title('num_snps_total vs num_coding_snps')
        plt.xlabel('num_snps_total')
        plt.ylabel('num_coding_snps')
        plt.savefig(f'{results_directory}/coding_snps_vs_total.png')
        plt.close()

        plt.figure(figsize=(8, 6))
        plt.scatter(non_coding_comparison['num_snps_total'], non_coding_comparison['num_non_coding'], color='b', alpha=0.6)
        plt.title('num_snps_total vs num_non_coding_snps')
        plt.xlabel('num_snps_total')
        plt.ylabel('num_non_coding_snps')
        plt.savefig(f'{results_directory}/non_coding_snps_vs_total.png')
        plt.close()

    return 

if __name__ == "__main__":
    main()
