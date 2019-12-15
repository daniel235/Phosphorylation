import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle

#import dev files
import alias

#read unique kinases
filename = "./data/pickles/uniqueKinases"
with open(filename, 'rb+') as f:
    uniqueKinases = pickle.load(f)


#family lookup
family_data = pd.read_csv("./data/kinaseClass.txt", delimiter=",")
#create list to get names of families in our phosphorylation data
family = []
for i in range(len(family_data)):
    #get unique family names
    if family_data['Classification'][i] not in family:
        family.append(family_data['Classification'][i])

    if len(family) > 11:
        break


kinase_dict = {}
alias_object = alias.Alias("./data/info_table.csv")
for kinase in uniqueKinases:
    #?add family 
    if kinase not in family_data['Gene'].tolist():
        #!get alias and try
        newKinase = alias_object.get_main_kinase(kinase)
        kinase_dict[kinase] = family_data['Classification'][family_data['Gene'].tolist().index(newKinase)]

    else:
        kinase_dict[kinase] = family_data['Classification'][family_data['Gene'].tolist().index(kinase)]
        


#create matrix
interaction_matrix = np.zeros(shape=(12, 12))


#get correlation values
filename = "./data/pickles/clusterGroups"
with open(filename, 'rb+') as f:
    cluster_groups = pickle.load(f)



#go through cluster groups
for typeNum in range(len(cluster_groups)):
    #required data structures 
    previousFams = []
    famDict = {}
    interaction_matrix = np.zeros(shape=(12,12))
    for j in range(len(cluster_groups[typeNum])):
        #inside the cluster group k
        clusterLen = len(cluster_groups[typeNum][j])
        print("cluster_groups ", cluster_groups[typeNum][j]) 
        for k in cluster_groups[typeNum][j]:
            #add instance to interaction matrix
            #get index from unique kinases
            #search for family kinase
            try:
                fam = family_data['Gene'].tolist().index(k) 
                family_found = True

            except ValueError:
                try:
                    fam = family_data['Gene'].tolist().index(alias_object.get_main_kinase(k))
                    family_found = True
                except:
                    family_found = False
                
            if(family_found):
                #set family index
                index = family.index(family_data['Classification'][fam])

                #?increment family in dictionary if it appears in cluster
                if family[index] in famDict:
                    famDict[family[index]] += 1
                else:
                    famDict[family[index]] = 1

        #!before moving to new cluster
        #iterate through famdict to start filling in interaction matrix
        for fam, cnt in famDict.items():
            for fam2, cnt2 in famDict.items():
                if cnt > (clusterLen * .2) and cnt2 > (clusterLen * .2) and fam != fam2:
                    interaction_matrix[family.index(fam), family.index(fam2)] += 1

        #reset famdict
        previousFams.append(famDict)
        famDict = {}

    #visualize data 
    fig, ax = plt.subplots()
    im = ax.imshow(interaction_matrix)

    ax.set_xticks(np.arange(12))
    ax.set_yticks(np.arange(12))
    ax.set_xticklabels(family)
    ax.set_yticklabels(family)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    for i in range(12):
        for j in range(12):
            text = ax.text(j, i, interaction_matrix[i, j], ha="center", va="center", color="w")

    ax.set_title("Interaction matrix")
    fig.tight_layout()
    plt.savefig(("./results/" + str(typeNum) + "InteractionMatrix.png"))
    plt.show()


#?class used for random interaction matrix
#?random interaction matrix
class InteractionMatrix:
    def __init__(self, clusterGroup, name=None):
        self.cluster_groups = clusterGroup
        self.kinases = None
        self.reps = 10
        self.family = None
        self.kinase_dict = None
        self.get_data()


    def run_interaction(self):
        kinase_dict = self.kinase_dict
        alias_object = alias.Alias("./data/info_table.csv")
        family = self.family
        #family lookup
        family_data = pd.read_csv("./data/kinaseClass.txt", delimiter=",")

        #create matrix
        interaction_matrix = np.zeros(shape=(12, 12))
        cluster_groups = self.cluster_groups
        
        #go through cluster groups
        #required data structures 
        famDict = {}
        for j in range(len(cluster_groups[0])):
            #inside the cluster group k
            clusterLen = len(cluster_groups[0][j])
            for k in cluster_groups[0][j]:
                #add instance to interaction matrix
                #get index from unique kinases
                #search for family kinase
                try:
                    fam = family_data['Gene'].tolist().index(k) 
                except ValueError:
                    fam = family_data['Gene'].tolist().index(alias_object.get_main_kinase(k))

                #set family index
                index = family.index(family_data['Classification'][fam])
            
                #?increment family in dictionary if it appears in cluster
                if family[index] in famDict:
                    famDict[family[index]] += 1
                else:
                    famDict[family[index]] = 1


            #!before moving to new cluster
            #iterate through famdict to start filling in interaction matrix
            for fam, cnt in famDict.items():
                for fam2, cnt2 in famDict.items():
                    print(fam, fam2)
                    if cnt > (clusterLen * .3) and cnt2 > (clusterLen * .3) and fam != fam2:
                        interaction_matrix[family.index(fam), family.index(fam2)] += 1

            famDict = {}

            
        self.save_matrix(interaction_matrix, family)
        
        return interaction_matrix


    def get_data(self):
        #read unique kinases
        filename = "./data/pickles/uniqueKinases"
        with open(filename, 'rb+') as f:
            uniqueKinases = pickle.load(f)


        #family lookup
        family_data = pd.read_csv("./data/kinaseClass.txt", delimiter=",")
        family = []
        for i in range(len(family_data)):
            #get unique family names
            if family_data['Classification'][i] not in family:
                family.append(family_data['Classification'][i])


        kinase_dict = {}
        alias_object = alias.Alias("./data/info_table.csv")
        for kinase in uniqueKinases:
            #?add family 
            if kinase not in family_data['Gene'].tolist():
                #!get alias and try
                newKinase = alias_object.get_main_kinase(kinase)
                kinase_dict[kinase] = family_data['Classification'][family_data['Gene'].tolist().index(newKinase)]

            else:
                kinase_dict[kinase] = family_data['Classification'][family_data['Gene'].tolist().index(kinase)]
                #check for family

        self.kinase_dict = kinase_dict
        self.family = family
        

    def save_matrix(self, interactionMatrixFig, family):
        #visualize data 
        fig, ax = plt.subplots()
        im = ax.imshow(interactionMatrixFig)

        ax.set_xticks(np.arange(12))
        ax.set_yticks(np.arange(12))
        ax.set_xticklabels(family)
        ax.set_yticklabels(family)

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        for i in range(12):
            for j in range(12):
                text = ax.text(j, i, interactionMatrixFig[i, j], ha="center", va="center", color="w")

        ax.set_title("Interaction matrix")
        fig.tight_layout()
        plt.savefig(("./results/RandomInteractionMatrix.png"))
        #plt.show()
        plt.close()
        return




