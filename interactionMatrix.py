import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle

import alias

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
    kinase_dict[kinase] = alias_object.get_main_kinase(kinase)
    #check for family



#create matrix
interaction_matrix = np.zeros(shape=(12, 12))

print(interaction_matrix)

#get correlation values
filename = "./data/pickles/clusterGroups"
with open(filename, 'rb+') as f:
    cg = pickle.load(f)


#go through cluster groups
for j in range(len(cg[0])):
    #inside the cluster group k
    for k in cg[0][j]:
        for l in cg[0][j]:
            if k != l:
                #add instance to interaction matrix
                #get index from unique kinases
                aliasKinase = kinase_dict[k]
                aliasKinase2 = kinase_dict[l]
                print(aliasKinase)
                print(aliasKinase2)
                #search for family kinase
                fam1 = family_data['Gene'].tolist().index(aliasKinase)
                fam2 = family_data['Gene'].tolist().index(aliasKinase2)
                index1 = family.index(family_data['Classification'][fam1])
                index2 = family.index(family_data['Classification'][fam2])
                interaction_matrix[index1, index2] += 1
                interaction_matrix[index2, index1] += 1


print(interaction_matrix)

#visualize data 
fig, ax = plt.subplots()
im = ax.imshow(interaction_matrix)

ax.set_xticks(np.arange(len(uniqueKinases)))
ax.set_yticks(np.arange(len(uniqueKinases)))
ax.set_xticklabels(uniqueKinases)
ax.set_yticklabels(uniqueKinases)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

for i in range(len(uniqueKinases)):
    for j in range(len(uniqueKinases)):
        text = ax.text(j, i, interaction_matrix[i, j], ha="center", va="center", color="w")

ax.set_title("Interaction matrix")
fig.tight_layout()
plt.show()

