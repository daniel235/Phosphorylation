import pandas as pd
import os
from pathlib import Path
import pickle

#change directory 
path = Path(os.getcwd())

#script to double check all kinases 
print("path", path)
try:
    data = pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object)
except: 
    path = path.parent
    os.chdir(str(path.resolve()))
    data = pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object)


#combine columns
for i in range(len(data['geneSymbol'])):
    data['geneSymbol'][i] = str(data['geneSymbol'][i]) + "-" + str(data['variableSites'][i])[:-2]


#run through all kinases
kinase_data = pd.read_csv("./data/KSA_human.txt", sep="\t")

matched_data = {}

for i in range(len(kinase_data)):
    kinase_data['Substrate'][i] = str(kinase_data['Substrate'][i]) + "-" + str(kinase_data['Site'][i])
    #look for substrate in phos file
    for j in range(len(data)):
        if kinase_data['Substrate'][i] == data['geneSymbol'][j]:
            if kinase_data['Kinase'][i] in matched_data:
                matched_data[kinase_data['Kinase'][i]].append(kinase_data['Substrate'][i])
            else:
                matched_data[kinase_data['Kinase'][i]] = [kinase_data['Substrate'][i]]
    
#pickle object
with open("./data/pickles/matchData", 'wb+') as f:
    pickle.dump(matched_data, f)


print(matched_data)

