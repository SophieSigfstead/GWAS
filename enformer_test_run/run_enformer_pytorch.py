from enformer_pytorch import from_pretrained, seq_indices_to_one_hot
import torch

def main():

    model = from_pretrained('EleutherAI/enformer-official-rough')

    seq = torch.randint(0, 5, (1, 196_608))

    one_hot = seq_indices_to_one_hot(seq)

    pred = model(one_hot)

    print(pred)

    return

if __name__ == "__main__":
    main()