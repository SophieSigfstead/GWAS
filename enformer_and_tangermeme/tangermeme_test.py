from enformer_pytorch import from_pretrained, seq_indices_to_one_hot
import torch
import tangermeme
from tangermeme.deep_lift_shap import deep_lift_shap
from matplotlib import pyplot as plt
from tangermeme.plot import plot_logo
from tangermeme.utils import random_one_hot
from tangermeme.ersatz import substitute
import sys
from matplotlib import pyplot as plt
import seaborn; seaborn.set_style('whitegrid')
from tangermeme.plot import plot_logo
from tangermeme.ersatz import substitute

class ModelWrapper:
    def __init__(self, model, head_key='human'):
        self.model = model
        self.head_key = head_key

    def __call__(self, *args, **kwargs):
        transposed_input = args[0].transpose(1, 2)
        preds = self.model(transposed_input)
        # extract the 'human' predictions 
        human_preds = preds[self.head_key]

        #human_preds_select_track = human_preds[:, :, 200]

        return human_preds
    
    def to(self, device):
        self.model = self.model.to(device)
        return self

    def eval(self):
        self.model.eval()
        return self

    def modules(self):
        return self.model.modules()
    
    def apply(self, fn):
        self.model.apply(fn)
        return self

def main():

    model = from_pretrained('EleutherAI/enformer-official-rough', use_tf_gamma = False)
    model = model.to(device = 'cpu')

    #one_hot = random_one_hot((1, 196608, 4)).type(torch.float32)
    one_hot = torch.randint(0, 4, (1, 196608))
    one_hot = seq_indices_to_one_hot(one_hot)


    wrapped_model = ModelWrapper(model)
    wrapped_pred = wrapped_model(one_hot.transpose(2,1))
    print(wrapped_pred.shape)
    del wrapped_pred

    X = random_one_hot((1, 4, 196608)).type(torch.float32)
    print(X.shape)
    reference_seq = random_one_hot((1, 4, 196608)).type(torch.float32)
    reference_seq = torch.unsqueeze(reference_seq, dim = 1)
    print(reference_seq.shape)
    wrapped_model_X = wrapped_model(X)
    print(wrapped_model_X.shape)
    del wrapped_model_X

    X = substitute(X, "GTGACTCATC")

    # empty the cache
    torch.cuda.empty_cache()

    #with torch.cuda.amp.autocast(): # this is supposed to reduce
    X_attr = deep_lift_shap(wrapped_model, X, target=200, device = 'cpu', random_state=0, batch_size = 16, references = reference_seq, verbose = True)

    plt.figure(figsize=(10, 2))
    ax = plt.subplot(111)
    plot_logo(X_attr[0, :, 950:1050], ax=ax)

    plt.xlabel("Genomic Coordinate")
    plt.ylabel("Attributions")
    plt.title("DeepLIFT Attributions for GM12878 ETS1")
    plt.ylim(-0.6, 0.6)
    plt.show()
    plt.savefig("image.png")


    return



if __name__ == "__main__":
    main()
    sys.exit()