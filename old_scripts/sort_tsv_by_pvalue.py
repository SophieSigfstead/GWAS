import sys
import pandas as pd
import os

def main(path_to_tsv, output_dir):

    # sort the tsv rows by the column p_value, from smallest to largest
    df = pd.read_csv(path_to_tsv, sep='\t')
    df = df.sort_values(by='p_value', ascending=True)
    
    # add a column for the ranking of the row, titled "ranking_of_p_value"
    df['ranking_of_p_value'] = range(1, len(df) + 1)
    
    # create the output file path
    base_name = os.path.basename(path_to_tsv)
    base_name_without_ext = os.path.splitext(base_name)[0]
    output_file = f"{base_name_without_ext}_sorted.tsv"
    output_path = os.path.join(output_dir, output_file)
    
    df.to_csv(output_path, sep='\t', index=False)
    print(f"File saved to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <path to tsv> <path to output dir>")
    else:
        path_to_tsv = sys.argv[1]
        output_dir = sys.argv[2]
        main(path_to_tsv, output_dir)
    sys.exit()

