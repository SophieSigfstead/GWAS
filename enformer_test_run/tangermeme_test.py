from enformer_pytorch import from_pretrained, seq_indices_to_one_hot
import torch
import tangermeme
from tangermeme.deep_lift_shap import deep_lift_shap
from matplotlib import pyplot as plt
from tangermeme.plot import plot_logo
from tangermeme.utils import random_one_hot
def main():

    model = from_pretrained('EleutherAI/enformer-official-rough')

    one_hot = random_one_hot((1, 4, 393216)).type(torch.float32)

    pred = model(one_hot)

    print(one_hot.shape)

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