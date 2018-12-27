import numpy as np
import pandas as pd
import tensor


class Pipe_line:
    def find_matching_data(self):
        #create data object for each line
        pass


    def get_data(self):
        kinaseData = pd.read_csv('./data/Kinase_Substrates.txt', delimiter="\t")

        proteinInteraction = pd.read_csv('./data/Protein_Protein_Interaction.txt', delimiter="\t")

        phosphorylation = pd.read_csv('./data/phosphorylation_data.txt', delimiter="\t")
        proteinExpression = pd.read_csv('./data/ProteinExpression_data.txt' ,delimiter="\t")
        #Phosphorylation name
        phosphoType = phosphorylation["Phosphosite"]

        proteinExpression = np.array(proteinExpression)
        phosphorylation = np.array(phosphorylation)
        proteinInteraction = np.array(proteinInteraction)

        print("prot size ", proteinInteraction.shape)
        print("phos size ", phosphorylation.shape)
        print("prot expre size ", proteinExpression.shape)

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
        key = ""
        keyTwo = ""
        for i in range(len(phosDataX) - 1):
            key = phosDataX[i][0]
            keyTwo = phosDataX[i + 1][0]
            if key != keyTwo:
                phosClass.append(keyTwo)

        for i in range(len(protExpressX) - 1):
            key = protExpressX[i][0]
            keyTwo = protExpressX[i + 1][0]
            if key != keyTwo:
                protClass.append(keyTwo)

        return [phosDataX, phosDataY, protExpressX, protExpressY, phosClass, protClass]