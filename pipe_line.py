import protein as pr
import numpy as np
import pandas as pd
import network


class Pipe_line:
    def __init__(self):
        self.kinaseData = pd.read_csv('./data/Kinase_Substrates.txt', delimiter="\t")
        self.proteinInteraction = pd.read_csv('./data/Protein_Protein_Interaction.txt', delimiter="\t")
        self.phosphorylation = pd.read_csv('./data/phosphorylation_data.txt', delimiter="\t")
        self.proteinExpression = pd.read_csv('./data/ProteinExpression_data.txt' ,delimiter="\t")
        

    def find_matching_data(self, data):
        proteinInteraction = np.array(self.proteinInteraction)

        #todo create data object for each line
        protein_names = []
        protein_objects = []
        id = 0
        p = None

        #check if protein is in phosphorylation/gene expression/and is a kinase substrate
        #first count check
        for i in range(len(data[0])):
            if data[0][i][0] not in protein_names:
                protein_names.append(data[0][i][0])
                p = pr.Protein(id)
                p.name_protein(data[0][i][0])
                p.count += 1
                protein_objects.append(p)
                id += 1

        #contains protein name in string
        substrate = self.kinaseData["Substrate"]
        substrate = np.array(substrate)
        newStr = ""
        #strip substrate sites
        for i in range(len(substrate)):
            for j in range(1, len(substrate[i])):
                if substrate[i][j] != '-':
                    newStr += substrate[i][j]
                else:
                    substrate[i] = str(newStr)
                    newStr = ""
                    break

        #turn np array into array to check for names
        sub_names = []
        for i in range(len(substrate)):
            sub_names.append("'" + substrate[i] + "'")


        #second count check
        for i in range(len(protein_names)):
            if protein_names[i] in sub_names:
                protein_objects[i].count += 1


        interaction_one = []
        interaction_two = []

        for i in range(len(proteinInteraction)):
            interaction_one.append(proteinInteraction[i][0])
            interaction_two.append(proteinInteraction[i][1])

        for i in range(len(interaction_one)):
            interaction_one[i] = "'" + str(interaction_one[i]) + "'"
            interaction_two[i] = "'" + str(interaction_two[i]) + "'"

        #third count check
        for i in range(len(protein_names)):
            if protein_names[i] in interaction_one:
                protein_objects[i].count += 1

            elif protein_names[i] in interaction_two:
                protein_objects[i].count += 1

        #data base are the protein objects that are in all 3 databases
        data_base = []
        for i in range(len(protein_names)):
            if protein_objects[i].count == 3:
                data_base.append(protein_objects[i])

        return data_base

    #function to retrieve data to be used in other files
    def get_data(self):
        #convert data into workable numpy arrays
        proteinExpression = np.array(self.proteinExpression)
        phosphorylation = np.array(self.phosphorylation)

        #phoshporylation data 
        phosDataX = []
        phosDataY = []

        #protein Expression 
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

        return [phosDataX, phosDataY, protExpressX, protExpressY, phosClass, protClass, self.kinaseData]


    #pass in kinase argument and return it's substrates from the 
    #phosphorylation data and expression data
    def grab_substrates(self, kinase):
        #load kinase indexes

        #transform kinase to string to match database kinases
        kinase = "'" + kinase + "'"
        kinaseColumn = self.kinaseData["Kinase"]
        list_of_substrates = []
        for i in range(len(kinaseColumn)): 
            if(kinaseColumn[i] == kinase):
                #add to substrate array
                list_of_substrates.append(self.kinaseData["Substrate"][i])

        return list_of_substrates
    
    def find_kinase(self, psite):
        k = self.kinaseData["Substrate"].tolist()
        try:
            ksubIndex = k.index("'" + psite + "'")
        except:
            ksubIndex = False

        return ksubIndex

    def strip_sites(self, site):
        # strip substrate sites
        #pass in single site
        newStr = ""
        index = 0
        for i in range(1, len(site)):
            #grab last hyphen
            if site[i] != '-':
                newStr += site[i]
            else:
                newStr += site[i]
                index = i

        #replace string up to index
        newStr = newStr[:index-1]


        return newStr
