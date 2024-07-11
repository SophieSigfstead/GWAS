import sys
import os
import h5py
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import gc

def process_file_in_chunks(file_path, combined_output_path, filtered_indices, chunk_size=10000):
    with h5py.File(file_path, 'r') as h5file:
        alt = h5file['alt']
        chr = h5file['chr']
        pos = h5file['pos']
        ref = h5file['ref']
        snp = h5file['snp']
        SAD = h5file['SAD']
        
        num_rows = alt.shape[0]

        # create an empty CSV for combined data
        headers = ['alt', 'chr', 'pos', 'ref', 'snp'] + [f'SAD{index}' for index in range(len(filtered_indices))]
        pd.DataFrame(columns=headers).to_csv(combined_output_path, index=False)

        for start in range(0, num_rows, chunk_size):
            end = min(start + chunk_size, num_rows)
            df_single_chunk = pd.DataFrame({
                'alt': alt[start:end],
                'chr': chr[start:end],
                'pos': pos[start:end],
                'ref': ref[start:end],
                'snp': snp[start:end]
            })

            df_SAD_chunk = pd.DataFrame(SAD[start:end, filtered_indices], columns=[f'SAD{index}' for index in range(len(filtered_indices))])

            df_combined_chunk = pd.concat([df_single_chunk, df_SAD_chunk], axis=1)

            # Append chunk to CSV
            df_combined_chunk.to_csv(combined_output_path, mode='a', header=False, index=False)

            print(f"Processed chunk {start} to {end} of {num_rows}")

            #gb collect
            gc.collect()

    return file_path

def main(input_path, output_path, chunk_size=10000):
    if not os.path.isdir(input_path):
        print(f"The input folder '{input_path}' does not exist.")
        return

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    file_names = [file_name for file_name in os.listdir(input_path) if file_name.endswith('.h5')]

    target_ids = None
    target_labels = None

    # Extract target data once and save
    for file_name in file_names:
        file_path = os.path.join(input_path, file_name)
        with h5py.File(file_path, 'r') as h5file:
            if target_ids is None:
                target_ids = h5file['target_ids'][:]
                target_labels = [label.decode('utf-8') for label in h5file['target_labels'][:]]
                # Save target data with numeric indices
                df_targets = pd.DataFrame({
                    'index': [i  for i in range(len(target_ids))],
                    'target_ids': target_ids,
                    'target_labels': target_labels
                })
                target_output_path = os.path.join(output_path, "targets.csv")
                df_targets.to_csv(target_output_path, index=False)
                print(f"Saved target data to {target_output_path}")
            break  # Only process the first file for target data

    # filter the indices of SAD columns to be saved
    filtered_indices = [i for i, label in enumerate(target_labels) if "DNASE" in label or "ATAC" in label]

    with ProcessPoolExecutor() as executor:
        futures = []
        for file_name in file_names:
            file_path = os.path.join(input_path, file_name)
            combined_output_path = os.path.join(output_path, f"{os.path.splitext(file_name)[0]}_combined.csv")
            
            futures.append(executor.submit(process_file_in_chunks, file_path, combined_output_path, filtered_indices, chunk_size))

        for future in futures:
            processed_file = future.result()
            print(f"Processed {processed_file}")

    print(f"All files processed. DataFrames saved in '{output_path}'.")
    return

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_path> <output_path>")
    else:
        input_path = sys.argv[1]
        output_path = sys.argv[2]
        main(input_path, output_path)
