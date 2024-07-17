import pandas as  pd
import matplotlib.pyplot as plt
import sys

def main(path_to_combined_csv):

    # the goal of this code is to establish a threshold for SNP inclusion 

    df = pd.read_csv(path_to_combined_csv)

    # count up how many have value for beta but not for SAD, and how many have value for SAD1 but not beta

    beta_not_SAD_count = df[(df['beta'].notna()) & (df['SAD1'].isna())].shape[0]
    SAD_not_beta_count = df[(df['SAD1'].notna()) & (df['beta'].isna())].shape[0]
    length_of_df = df.shape[0]

    print(f"Rows from tsv with no match in 1000 genomes: {beta_not_SAD_count}")
    print(f"Rows with from 1000 genomes with no match in tsv: {SAD_not_beta_count}")
    print(f"Length of df: {length_of_df}")
    
    # generate a plot of abs_average_SAD distribution of all rows
    plt.figure(figsize=(10, 6))
    plt.xlim(0, 0.0005)
    df['abs_average_SAD'].dropna().plot(kind='hist', bins=2000, edgecolor='black')
    plt.title('Distribution of abs_average_SAD')
    plt.xlabel('abs_average_SAD')
    plt.ylabel('Frequency')
    plt.savefig('./SNP_threshold_data_exploration/abs_average_SAD_distribution.png')

    plt.figure(figsize=(10, 6))
    plt.xlim(-0.0005, 0.0005)
    df['SAD1'].dropna().plot(kind='hist', bins=2000, edgecolor='black')
    plt.title('Distribution of SAD1')
    plt.xlabel('SAD1')
    plt.ylabel('Frequency')
    plt.savefig('./SNP_threshold_data_exploration/SAD1_distribution.png')

    plt.figure(figsize=(10, 6))
    plt.xlim(-0.0005, 0.0005)
    df['SAD9'].dropna().plot(kind='hist', bins=2000, edgecolor='black')
    plt.title('Distribution of SAD9')
    plt.xlabel('SAD9')
    plt.ylabel('Frequency')
    plt.savefig('./SNP_threshold_data_exploration/SAD9_distribution.png')

    # filter the combined csv file to include only rows that have values for beta and average_SAD

    filtered_df = df[df['beta'].notna() & df['average_SAD'].notna()]
    print(len(filtered_df))
    # generate a plot of correlation between beta and average_SAD

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['beta'], filtered_df['abs_average_SAD'], s = 0.25)
    plt.title('Correlation between beta and average_SAD')
    plt.xlabel('beta')
    plt.ylabel('average_SAD')
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between beta and abs_average_SAD.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['beta'], filtered_df['SAD1'], s = 0.25)
    plt.title('Correlation between beta and SAD1')
    plt.xlabel('beta')
    plt.ylabel('SAD1')
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between beta and SAD1.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['beta'], filtered_df['SAD9'], s = 0.25)
    plt.title('Correlation between beta and SAD9')
    plt.xlabel('beta')
    plt.ylabel('SAD9')
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between beta and SAD9.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['p_value'], filtered_df['average_SAD'], s = 0.25)
    plt.title('Correlation between beta and average_SAD')
    plt.xlabel('p_value')
    plt.ylabel('average_SAD')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between pval and average_SAD.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['p_value'], filtered_df['SAD1'], s = 0.25)
    plt.title('Correlation between beta and SAD1')
    plt.xlabel('p_value')
    plt.ylabel('SAD1')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between pval and SAD1.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['p_value'], filtered_df['SAD9'], s = 0.25)
    plt.title('Correlation between beta and SAD9')
    plt.xlabel('p_value')
    plt.ylabel('SAD9')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between pval and SAD9.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df['p_value'], filtered_df['abs_average_SAD'], s = 0.25)
    plt.title('Correlation between beta and SAD9')
    plt.xlabel('p_value')
    plt.ylabel('average_SAD')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/Correlation between pval and abs_average_SAD.png')


    df_pval = pd.read_csv('./sorted_GWAS_summary_data/GCST90277450_sorted.tsv', sep='\t')
    filtered_df_pval = df_pval[df_pval['p_value'] < 5e-8]
    filtered_file_path = './sorted_GWAS_summary_data/GCST90277450_sorted_significant.csv'
    filtered_df_pval.to_csv(filtered_file_path, sep='\t', index=False)

    filtered_df_significant = df[df['beta'].notna() & df['average_SAD'].notna() & (df['p_value'] < 5e-8) & df['p_value'].notna()]
    print('length filtered df significant')
    print(len(filtered_df_significant))

    # generate a plot of abs_average_SAD distribution of all rows
    plt.figure(figsize=(10, 6))
    plt.xlim(0, 0.0005)
    filtered_df_significant['abs_average_SAD'].dropna().plot(kind='hist', bins=2000, edgecolor='black')
    plt.title('Distribution of abs_average_SAD')
    plt.xlabel('abs_average_SAD')
    plt.ylabel('Frequency')
    plt.savefig('./SNP_threshold_data_exploration/filtered abs_average_SAD_distribution.png')

    plt.figure(figsize=(10, 6))
    plt.xlim(-0.0005, 0.0005)
    filtered_df_significant['SAD1'].dropna().plot(kind='hist', bins=2000, edgecolor='black')
    plt.title('Distribution of SAD1')
    plt.xlabel('SAD1')
    plt.ylabel('Frequency')
    plt.savefig('./SNP_threshold_data_exploration/filtered SAD1_distribution.png')

    plt.figure(figsize=(10, 6))
    plt.xlim(-0.0005, 0.0005)
    filtered_df_significant['SAD9'].dropna().plot(kind='hist', bins=2000, edgecolor='black')
    plt.title('Distribution of SAD9')
    plt.xlabel('SAD9')
    plt.ylabel('Frequency')
    plt.savefig('./SNP_threshold_data_exploration/filtered SAD9_distribution.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['beta'], filtered_df_significant['abs_average_SAD'], s = 0.25)
    plt.title('Correlation between beta and average_SAD')
    plt.xlabel('beta')
    plt.ylabel('average_SAD')
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between beta and abs_average_SAD.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['beta'], filtered_df_significant['SAD1'], s = 0.25)
    plt.title('Correlation between beta and SAD1')
    plt.xlabel('beta')
    plt.ylabel('SAD1')
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between beta and SAD1.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['beta'], filtered_df_significant['SAD9'], s = 0.25)
    plt.title('Correlation between beta and SAD9')
    plt.xlabel('beta')
    plt.ylabel('SAD9')
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between beta and SAD9.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['p_value'], filtered_df_significant['average_SAD'], s = 0.25)
    plt.title('Correlation between beta and average_SAD')
    plt.xlabel('p_value')
    plt.ylabel('average_SAD')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between pval and average_SAD.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['p_value'], filtered_df_significant['SAD1'], s = 0.25)
    plt.title('Correlation between beta and SAD1')
    plt.xlabel('p_value')
    plt.ylabel('SAD1')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between pval and SAD1.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['p_value'], filtered_df_significant['SAD9'], s = 0.25)
    plt.title('Correlation between beta and SAD9')
    plt.xlabel('p_value')
    plt.ylabel('SAD9')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between pval and SAD9.png')

    plt.figure(figsize=(10, 6))
    plt.scatter(filtered_df_significant['p_value'], filtered_df_significant['abs_average_SAD'], s = 0.25)
    plt.title('Correlation between beta and SAD9')
    plt.xlabel('p_value')
    plt.ylabel('average_SAD')
    plt.xlim(0, 5e-8)
    plt.ylim(-0.0005, 0.0005)
    plt.grid(True)
    plt.savefig('./SNP_threshold_data_exploration/filtered Correlation between pval and abs_average_SAD.png')









    return


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path for combined csv>")
    else:
        path_to_combined_csv = sys.argv[1]
        main(path_to_combined_csv)
    sys.exit()