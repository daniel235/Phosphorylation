import pca
import numpy as np
import pandas as pd
import clean


phos_data = pd.read_csv("./data/phosphorylation_data.txt", delim_whitespace=True)
phos_data = np.array(phos_data)
for i in range(len(phos_data)):
    phos_data[i][0] = phos_data[i][0].replace("'", '')

sub_data = pd.read_csv("./data/KSA_human.txt", delim_whitespace=True)
sub_data = np.array(sub_data)
for i in range(len(sub_data)):
    sub_data[i][1] = sub_data[i][1] + "-" + str(sub_data[i][2])
    

sub_data = np.delete(sub_data, 2, 1)


set1 = set(list(phos_data[:,0]))
set2 = set(list(sub_data[:,1]))


print(set1.intersection(set2))