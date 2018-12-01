import tensorflow as tf
import numpy as np
import pandas as pd


kinaseData = pd.read_csv('./data/Kinase_Substrates.txt', delimiter="\t")

phosphorylation = pd.read_csv('./data/phosphorylation_data.txt', delimiter="\t")

#Phosphorylation name
phosphoType = phosphorylation["Phosphosite"]

phosphorylation = np.array(phosphorylation)

phosDataX = []
phosDataY = []


#phosphorylation data separated with xname and x data
for i in range(len(phosphorylation)):
    for j in range(1, 19):
        phosDataX.append([phosphorylation[i][0], phosphorylation[i][j]])


#need to link y data








