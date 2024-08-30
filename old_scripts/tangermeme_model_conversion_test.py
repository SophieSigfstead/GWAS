import tensorflow as tf
import sys
import numpy as np
import tensorflow_hub as hub

def save_tf_model():
    enformer_model = hub.load(
        "https://kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model
    enformer_model.save('saved_enformer_model_tf')
    return
def convert_to_onnx():
    return

def convert_to_pytorch():
    return

def verify_conversion():
    return


def main():
    save_tf_model()
    convert_to_onnx()
    convert_to_pytorch()
    verify_conversion()
    return


if __name__ == "__main__":
    main()
    sys.exit()