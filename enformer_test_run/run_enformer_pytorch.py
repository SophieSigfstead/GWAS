from enformer_pytorch import from_pretrained, seq_indices_to_one_hot
import torch
import sys

def main():

    model = from_pretrained('EleutherAI/enformer-official-rough')

    # creates sequence of numbers (0, 1, 2, 3)
    seq = torch.randint(0, 4, (1, 196608))
    print(seq.shape)

    one_hot = seq_indices_to_one_hot(seq)

    pred = model(one_hot)

    print(pred)

    return

if __name__ == "__main__":
    main()
    sys.exit()