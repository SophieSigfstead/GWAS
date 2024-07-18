import sys
import pandas as pd
import numpy

def threshold_calculator(df, sd, track):
    df_with_sad = df[df[track].notna()]
    print("len df with sad", len(df_with_sad))
    std_dev = df_with_sad[track].std().astype(float)
    print("std_dev", std_dev)
    mean = df_with_sad[track].mean().astype(float)
    print("mean:", mean)
    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)
    return upper_threshold, lower_threshold

def main(path_to_combined_csv, sd, comparison_study):

    TRACK_LIST = ['SAD1', 'SAD9']
    sd = float(sd)

    df = pd.read_csv(path_to_combined_csv)
    n_snps_prev = len(df)
    print('n_snps_prev:', n_snps_prev)

    for track in TRACK_LIST:
        df[track] = df[track].astype(float)
        print(df[track][0:5])
    
    
    filtered_df = df[df[TRACK_LIST[0]].notna()]

    output_columns = ['SAD1', 'SAD9', 'alt', 'chr', 'pos', 'ref', 'snp', 'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'average_SAD', 'abs_average_SAD', 'ranking']
    result_df = pd.DataFrame(columns=output_columns)

    threshold_dict = {}
    


    for track in TRACK_LIST:
        upper_threshold, lower_threshold = threshold_calculator(filtered_df, sd, track)
        threshold_dict[track] = {'upper': upper_threshold, 'lower': lower_threshold}

    print(threshold_dict)

    for track in TRACK_LIST:
        for chr_num in range(1, 23):
            chr_df = filtered_df[filtered_df['chr'] == chr_num]
            chr_df = chr_df[(chr_df['average_SAD'] > threshold_dict[track]['upper']) | (chr_df['average_SAD'] < threshold_dict[track]['lower'])]
            chr_df = chr_df.sort_values(by='p_value', ascending=True).reset_index(drop=True)

            while not chr_df.empty:
                top_snp = chr_df.iloc[0]
                result_df = pd.concat([result_df, pd.DataFrame([top_snp], columns=output_columns)], ignore_index=True)
                pos_sig_snp = top_snp['pos']
                chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + 1000000) & (chr_df['pos'] > pos_sig_snp - 1000000))]
                chr_df = chr_df.reset_index(drop=True)

    print("len result df", len(result_df))
    
    # modification: take the set of the result df rows
    result_df = pd.DataFrame(list(set([tuple(row) for row in result_df.to_records(index=False)])), columns=output_columns)

    n_snps_filtered = len(result_df)

    # threshold result according to new pvalues
    new_p_value_threshold = 5e-8 * (n_snps_prev / n_snps_filtered)
  
    print('new_p_value_threshold:', new_p_value_threshold)

    result_df = result_df[result_df['p_value'] < new_p_value_threshold]

    output_file = f'./filtered_snps_threshold={sd}.csv'
    result_df.to_csv(output_file, index=False)

    print(f"Length of the significant snps file: {len(result_df)}")

    # Compare with mdd_sig_snps.csv
    mdd_sig_snps = pd.read_csv('mdd_sig_snps.csv')
    matching_snps = pd.merge(result_df, mdd_sig_snps, left_on=['chr', 'pos'], right_on=['chromosome', 'base_pair_location'])
    print(f"Number of matching rows with mdd_sig_snps.csv: {len(matching_snps)}")
    return

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <path to combined csv> <number of sd for threshold> <comparison study>")
    else:
        path_to_combined_csv = sys.argv[1]
        sd = sys.argv[2]
        comparison_study = sys.argv[3]
        main(path_to_combined_csv, sd, comparison_study)
    sys.exit()

