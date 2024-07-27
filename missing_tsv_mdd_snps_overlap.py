import sys
import pandas as pd

def main(path_to_tsv,gwas_snps):


    tsv_df = pd.read_csv(path_to_tsv)
    tsv_df = tsv_df[tsv_df['SAD669'].isna()]
    print(len(tsv_df))
    pvalue_low = tsv_df[tsv_df['P'] < 5e-8]
    print(pvalue_low)
    print(pvalue_low.head())

    gwas_sig_snps = pd.read_csv(gwas_snps)
    print(gwas_sig_snps.head())

    matching_snps = pd.merge(tsv_df, gwas_sig_snps, left_on=['MarkerName'], right_on=['Marker_Name'])
    matching_snps.to_csv('./matching_snps_result')
    print(f"Number of matching rows with mdd_sig_snps.csv: {len(matching_snps)}")

    return


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <path for tsv missing rows from 1000 genomes> <path of gwas mdd sig snps>")
    else:
        path_to_tsv = sys.argv[1]
        gwas_snps = sys.argv[2]
        main(path_to_tsv, gwas_snps)
    sys.exit()