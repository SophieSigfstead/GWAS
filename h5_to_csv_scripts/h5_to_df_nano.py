import sys
import os
import h5py
import pandas as pd
import gc

def process_file_in_chunks(file_path, combined_output_path, target_output_path, chunk_size=10000):
    with h5py.File(file_path, 'r') as h5file:
        alt = h5file['alt']
        chr = h5file['chr']
        pos = h5file['pos']
        ref = h5file['ref']
        snp = h5file['snp']
        SAD = h5file['SAD']
        SAR = h5file['SAR']
        target_ids = h5file['target_ids'][:]
        target_labels = h5file['target_labels'][:]

        num_rows = alt.shape[0]

        # Create an empty CSV for combined data with headers
        pd.DataFrame(columns=[
            'alt', 'chr', 'pos', 'ref', 'snp',
            *[f'SAD{i + 1}' for i in range(SAD.shape[1])],
            *[f'SAR{i + 1}' for i in range(SAR.shape[1])]
        ]).to_csv(combined_output_path, index=False)

        for start in range(0, num_rows, chunk_size):
            end = min(start + chunk_size, num_rows)
            df_single_chunk = pd.DataFrame({
                'alt': alt[start:end],
                'chr': chr[start:end],
                'pos': pos[start:end],
                'ref': ref[start:end],
                'snp': snp[start:end]
            })

            df_SAD_chunk = pd.DataFrame(SAD[start:end], columns=[f'SAD{i + 1}' for i in range(SAD.shape[1])])
            df_SAR_chunk = pd.DataFrame(SAR[start:end], columns=[f'SAR{i + 1}' for i in range(SAR.shape[1])])

            df_combined_chunk = pd.concat([df_single_chunk, df_SAD_chunk, df_SAR_chunk], axis=1)

            # Append chunk to CSV
            df_combined_chunk.to_csv(combined_output_path, mode='a', header=False, index=False)

            # Print chunk processed message
            print(f"Processed chunk {start} to {end} of {num_rows}")

            # Explicitly run garbage collection
            gc.collect()

        # Save target data
        df_targets = pd.DataFrame({
            'target_ids': target_ids,
            'target_labels': target_labels
        })
        df_targets.to_csv(target_output_path, index=False)

        # Explicitly run garbage collection
        gc.collect()

    return file_path

def main(input_path, output_path, chunk_size=10000):
    if not os.path.isdir(input_path):
        print(f"The input folder '{input_path}' does not exist.")
        return

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    file_names = [file_name for file_name in os.listdir(input_path) if file_name.endswith('.h5')]

    for file_name in file_names:
        file_path = os.path.join(input_path, file_name)
        combined_output_path = os.path.join(output_path, f"{os.path.splitext(file_name)[0]}_combined.csv")
        target_output_path = os.path.join(output_path, f"{os.path.splitext(file_name)[0]}_targets.csv")

        processed_file = process_file_in_chunks(file_path, combined_output_path, target_output_path, chunk_size)
        print(f"Processed {processed_file}")

    print(f"All files processed. DataFrames saved in '{output_path}'.")
    return

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python main.py <input_folder> <output_folder> <chunk_size>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    chunk_size = int(sys.argv[3])
    main(input_path, output_path, chunk_size)
