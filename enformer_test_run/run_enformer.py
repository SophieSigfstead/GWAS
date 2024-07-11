import tensorflow as tf
import sys
import numpy as np
import tensorflow_hub as hub
def one_hot_encode_dna(sequence):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    one_hot_encoded = np.zeros((len(sequence), 4))
    for i, nucleotide in enumerate(sequence):
        one_hot_encoded[i, mapping[nucleotide]] = 1
    return one_hot_encoded

def main(model_path):

    SEQ_LENGTH = 393_216
    BATCH_SIZE = 1000

    enformer_model = hub.load("https://kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model

    dna_sequence = "AGCT" * (SEQ_LENGTH // 4)

    if len(dna_sequence) != SEQ_LENGTH:
        dna_sequence = dna_sequence[:SEQ_LENGTH]

    input_data = one_hot_encode_dna(dna_sequence)
    print("Input array shape: ", input_data.shape)

    input_array = np.array([input_data] * BATCH_SIZE, dtype=np.float32)
    print("Input array shape: ", input_array.shape)

    predictions = enformer_model.predict_on_batch(input_array)
    print(predictions)

    return

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_enformer.py <enformer_model_path>")
        sys.exit(1)

    model_path = sys.argv[1]
    main(model_path)


