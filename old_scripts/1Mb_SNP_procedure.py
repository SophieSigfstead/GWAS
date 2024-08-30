import pandas as pd
import sys

def threshold_calculator(df, sd):
    df_with_sad = df[df['SAD1'].notna()]
    std_dev = df_with_sad['average_SAD'].std()
    mean = df_with_sad['average_SAD'].mean()
    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)
    return upper_threshold, lower_threshold

def main(path_to_combined_csv, sd):

    sd = float(sd)

    df = pd.read_csv(path_to_combined_csv)

    UPPER_THRESHOLD, LOWER_THRESHOLD = threshold_calculator(df, sd)

    print("Number of SD for threshold:", sd)

    n_snps_prev = len(df)

    df = df[df['SAD1'].notna()]

    print('n_snps_prev:', n_snps_prev)

    output_columns = ['SAD1', 'SAD9', 'alt', 'chr', 'pos', 'ref', 'snp', 'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'average_SAD', 'abs_average_SAD', 'ranking']
    result_df = pd.DataFrame(columns=output_columns)

    n_snps_filtered = len(df[(df['average_SAD'] > UPPER_THRESHOLD) | (df['average_SAD'] < LOWER_THRESHOLD)])
    print('n_snps_filtered:', n_snps_filtered)

    for chr_num in range(1, 23):
        chr_df = df[df['chr'] == chr_num]
        chr_df = chr_df[(chr_df['average_SAD'] > UPPER_THRESHOLD) | (chr_df['average_SAD'] < LOWER_THRESHOLD)]
        chr_df = chr_df.sort_values(by='p_value', ascending=True).reset_index(drop=True)

        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = pd.concat([result_df, pd.DataFrame([top_snp], columns=output_columns)], ignore_index=True)
            pos_sig_snp = top_snp['pos']
            chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + 1000000) & (chr_df['pos'] > pos_sig_snp - 1000000))]
            chr_df = chr_df.reset_index(drop=True)

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
    if len(sys.argv) != 3:
        print("Usage: python script.py <path to combined csv> <number of sd for threshold>")
    else:
        path_to_combined_csv = sys.argv[1]
        sd = sys.argv[2]
        main(path_to_combined_csv, sd)
    sys.exit()