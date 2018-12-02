import tensorflow as tf
import numpy as np
import pandas as pd


kinaseData = pd.read_csv('./data/Kinase_Substrates.txt', delimiter="\t")

phosphorylation = pd.read_csv('./data/phosphorylation_data.txt', delimiter="\t")
proteinExpression = pd.read_csv('./data/ProteinExpression_data.txt' ,delimiter="\t")
#Phosphorylation name
phosphoType = phosphorylation["Phosphosite"]

proteinExpression = np.array(proteinExpression)
phosphorylation = np.array(phosphorylation)

phosDataX = []
phosDataY = []

protExpressX = []
protExpressY = []

phosClass = []
protClass = []
#phosphorylation data separated with xname and x data and added array of corresponding y values
for i in range(len(phosphorylation)):
    for j in range(1, 19):
        phosDataX.append([phosphorylation[i][0], phosphorylation[i][j]])
        if(i < 4037):
            protExpressX.append([proteinExpression[i][0], proteinExpression[i][j]])
            if j < 10:
                protExpressY.append("Basal")

            else:
                protExpressY.append("Luminal")
        if j < 10:
            phosDataY.append("Basal")
        else:
            phosDataY.append("Luminal")

#pass data to network
#prepare categories

for i in range(len(phosDataX)):
    if(i < len(phosDataX) - 1):
        print(phosDataX[i + 1][0], phosDataX[i][0])
        if phosDataX[i + 1][0] != phosDataX[i][0]:
            phosClass.append(phosDataX[i + 1][0])


print(len(phosClass))



