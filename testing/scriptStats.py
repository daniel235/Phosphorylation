import pandas as pd
import os
from pathlib import Path
import pickle

#change directory 
path = Path(os.getcwd())

'''script to double check if all kinases are showing up in phospho data'''
#*required data structures
matched_data = {}
greatest = 0


#?import phosphorylation data
#fixing path issue here
try:
    data = pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object)
#if not in correct path move to parent path to grab data
except: 
    path = path.parent
    os.chdir(str(path.resolve()))
    data = pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object)


#combine columns
for i in range(len(data['geneSymbol'])):
    data['geneSymbol'][i] = str(data['geneSymbol'][i]) + "-" + str(data['variableSites'][i])[:-2]


#?import kinase data
#run through all kinases
kinase_data = pd.read_csv("./data/KSA_human.txt", sep="\t")


#check if already have data for kinases
list_of_files = os.listdir("./data/pickles/")

#look for matchData files 
for i in range(len(list_of_files)):
    if "matchData" in list_of_files[i] and len(list_of_files[i]) > 10:
        #iterate to get the last file 
        if int(list_of_files[i][9:]) > greatest:
            greatest = int(list_of_files[i][9:])


for i in range(greatest + 1, len(kinase_data)):
    #every 3000 iterations pickle data
    if i % 3000 == 0 and i != 0:
        string = "./data/pickles/matchData" + str(i)
        with open(string, 'wb+') as f:
            pickle.dump(matched_data, f)


    kinase_data['Substrate'][i] = str(kinase_data['Substrate'][i]) + "-" + str(kinase_data['Site'][i])
    #look for substrate in phos file
    current_substrate = None
    for j in range(len(data)):
        if kinase_data['Substrate'][i] == data['geneSymbol'][j]:
            if kinase_data['Kinase'][i] in matched_data:
                matched_data[kinase_data['Kinase'][i]].append(kinase_data['Substrate'][i])
            else:
                matched_data[kinase_data['Kinase'][i]] = [kinase_data['Substrate'][i]]
    
            break

        


#pickle object
with open("./data/pickles/matchDataComplete", 'wb+') as f:
    pickle.dump(matched_data, f)



