from enformer_pytorch import from_pretrained, seq_indices_to_one_hot
import torch
import tangermeme
from tangermeme.deep_lift_shap import deep_lift_shap

def main():

    model = from_pretrained('EleutherAI/enformer-official-rough')

    seq = torch.randint(0, 5, (1, 196_608))

    one_hot = seq_indices_to_one_hot(seq)

    pred = model(one_hot)

    print(pred)

    X_attr = deep_lift_shap(model, seq, target='human', random_state=0)

    print(X_attr.shape)

    return


if __name__ == "__main__":
    main()
    sys.exit()