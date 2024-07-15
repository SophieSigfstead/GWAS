import os
import kagglehub
import kaggle
import shutil

path = kagglehub.model_download('deepmind/enformer/TensorFlow2/enformer/1')

print(path)

cwd =  os.getcwd()

file_name = os.path.basename(path)

destination = os.path.join(cwd, file_name)

shutil.move(path, destination)

#kagglehub.model_download('google/bert/tensorFlow2/answer-equivalence-bem')
