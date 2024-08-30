import sys
import pandas as pd
import os

def main(path_to_csv):

    # sort the tsv rows by the column p_value, from smallest to largest
    df = pd.read_csv(path_to_csv)
    df['SAD1abs'] = df['SAD1'].abs()
    df['SAD9abs'] = df['SAD9'].abs()
    SAD1df = df.sort_values(by='SAD1abs', ascending=False)
    SAD9df = df.sort_values(by='SAD9abs', ascending=False)
    
    # add a column for the ranking of the row, titled "ranking_of_p_value"
    SAD1df['ranking_of_SAD1'] = range(1, len(df) + 1)
    SAD9df['ranking_of_SAD9'] = range(1, len(df) + 1)
    
    SAD1df.to_csv('./SAD1sorted.csv',  index=False)

    SAD9df.to_csv('./SAD9sorted.csv',  index=False)

    return



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path to tsv> <path to output dir>")
    else:
        path_to_tsv = sys.argv[1]
        main(path_to_tsv)
    sys.exit()
