from enformer_pytorch import from_pretrained, seq_indices_to_one_hot
import torch
import tangermeme
from tangermeme.deep_lift_shap import deep_lift_shap
from matplotlib import pyplot as plt

def main():

    model = from_pretrained('EleutherAI/enformer-official-rough')

    seq = torch.randint(0, 5, (1, 196_608))

    one_hot = seq_indices_to_one_hot(seq)

    pred = model(one_hot)

    print(pred)

    X_attr = deep_lift_shap(model, one_hot, target='human', random_state=0)

    print(X_attr.shape)

    # Visualize attributions
    plt.figure(figsize=(10, 2))
    ax = plt.subplot(111)
    plot_logo(X_attr[0, :, 950:1050], ax=ax)

    plt.xlabel("Genomic Coordinate")
    plt.ylabel("Attributions")
    plt.title("DeepLIFT Attributions for Human")
    plt.ylim(-0.05, 0.35)
    plt.show()

    return


if __name__ == "__main__":
    main()
    sys.exit()