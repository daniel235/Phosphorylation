import pandas as pd
import numpy as np
import random
import clean
import cluster_data
import pca
import hierarchical
import math
import threading
import seaborn as sns
import interactionMatrix
import pickle
import compare_clusters
#import data

#file functions
import random_functions


#!todo fix scores bug
scores = []
#import data
cancer_data = np.array(pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object))

#clean data
cleanData = clean.cleanMatrix()
cleanData.data = cancer_data
cleanData.set_gene_site_column([0,1], True)
cleanData.omit_columns([16,17,18])
cleanData.clean_rows()


cancer_data = cleanData.data
im = np.zeros(shape=(12,12))

data_pipeline = cluster_data.PrepareClusterData("./data/KSA_human.txt")
data_pipeline.pfile = "./data/BreastCancerData.xlsx"

#fix kinases
#kinase alias name fix here
data_pipeline.convert_kinases("./data/KSA_human.txt")
tempKinases = np.array(pd.read_csv("./results/newPhosKinaseFile.txt"))
#fix kinase substrates columns
for i in range(len(data_pipeline.phosphositePlusKinaseData[:,1])-1):
    data_pipeline.phosphositePlusKinaseData[i,0] = tempKinases[i]
    data_pipeline.phosphositePlusKinaseData[i,1] = str(data_pipeline.phosphositePlusKinaseData[i,1]) + "-" + str(data_pipeline.phosphositePlusKinaseData[i,2])
    
data_pipeline.phosphositePlusKinaseData = data_pipeline.phosphositePlusKinaseData[:,0:-1]
#get unique kinases
data_pipeline.create_unique_kinases()

HierCluster = hierarchical.Hierarchical()

#add cluster to hyper geometric class
comparativeClusterGroups = compare_clusters.CompareCluster()
comparativeClusterGroups.setMainCluster()

#separate out psite column from numeric data to impute   
random_data = cancer_data[:,1:]

for indy in range(20):
    print("run ", indy)
    #impute random data
    #shuffle data in each row
    for i in range(len(random_data)):
        np.random.shuffle(random_data[i])
    #shuffle rows 
    np.random.shuffle(random_data)
    
    #put back imputed data into original cancer data
    for cell in range(len(cancer_data)):
        cancer_data[cell,1:] = random_data[cell]
        

    data_pipeline.CancerData = cancer_data
    data_pipeline.replace_with_average()

    #get substrate matrixes 
    subMatrix = data_pipeline.get_kinase_substrate_matrixes(2)

    #apply pca to my matrix
    X = []
    labels = []
    Xpoor = []
    labelsPoor = []
    Xrich = []
    labelsRich = []

    kinaseFeatures, poorKFeats, richKFeats, pfile = pca.getSVDdata("./data/KSA_human.txt", 10, obs="somethin", pfile=data_pipeline.fileName, matrix=subMatrix)
    for kinase, vector in kinaseFeatures.items():
        X.append(vector)
        labels.append(kinase)

    for kinase, vector in poorKFeats.items():
        Xpoor.append(vector)
        labelsPoor.append(kinase)

    for kinase, vector in richKFeats.items():
        Xrich.append(vector)
        labelsRich.append(kinase)
    
    #add data to hiercluster
    HierCluster.X = X
    HierCluster.labels = labels
    HierCluster.kinaseFeatures = kinaseFeatures
    HierCluster.poorKFeats = poorKFeats
    HierCluster.richKFeats = richKFeats
    HierCluster.labelsPoor = labelsPoor
    HierCluster.labelsRich = labelsRich
    HierCluster.pfile = pfile
    HierCluster.clusterMethod("notpca", pfile)

    #filter out kinases not in phosphorylation data
    comparativeClusterGroups.filter_phospho_kinases(labels)
    clLen = 0
    #get number of non empty family nodes
    family_indexes = []
    for i in range(len(comparativeClusterGroups.all_cluster_nodes[0])):
        if len(comparativeClusterGroups.all_cluster_nodes[0][i].data) != 0:
            family_indexes.append(i)
            clLen += 1

    HierCluster.start_hierarchical_clustering(clLen)
    comparativeClusterGroups.add_cluster(HierCluster.clusters)

    filename = "./data/pickles/randomClusterGroups"
    with open(filename, 'wb+') as f:
        pickle.dump(comparativeClusterGroups.all_clusters, f)

    comparativeClusterGroups.create_graph()
    comparativeClusterGroups.get_edge_scores(random=True)
    comparativeClusterGroups.display_stats()
    #todo get average of scores
    #!in range 12

    #?start random interaction
    matrix = interactionMatrix.InteractionMatrix(comparativeClusterGroups.all_clusters)
    im = np.add(im, matrix.run_interaction())
    
    counter = 0
    for ik in range(clLen): 
        for k in range(clLen): #12
            if indy == 0:
                if ik != k:
                    scores.append(comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[0][k].name])
                else:
                    scores.append(0)

            else:
                scores[counter] += comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[0][k].name]

            counter += 1

    
#divide scores by 50
for i in range(len(scores)):
    scores[i] = scores[i] / 20


#averaging interaction matrix
im = np.floor_divide(im, 20)
matrix.save_matrix(im, matrix.family)


##pickle data
filename="data/pickles/randomScores"
with open(filename, 'wb+') as f:
    pickle.dump(scores, f)

print("current scores ", scores)

with open("./data/pickles/randomScores", 'rb+') as f:
    obj = pickle.load(f)

with open("./data/pickles/clusterNodes", 'rb+') as f:
    cgNodes = pickle.load(f)


random_functions.get_sig_scores(obj, cgNodes)
#distro_scores(obj)

