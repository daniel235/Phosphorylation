import pandas as pd
import numpy as np
import random
import clean
import cluster_data
import pca
import hierarchical
import threading
import pickle
import compare_clusters
#import data

#!todo fix scores bug
scores = []


cancer_data = np.array(pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object))


#clean data
cleanData = clean.cleanMatrix()
cleanData.data = cancer_data
cleanData.set_gene_site_column([0,1], True)
cleanData.omit_columns([16,17,18])
cleanData.clean_rows()

cancer_data = cleanData.data


for indy in range(100):
    print("run ", indy)
    #insert random data into cells
    for i in range(len(cancer_data)):
        for j in range(1, len(cancer_data[i])):
            cancer_data[i][j] = random.uniform(-2, 3)
            #todo need to get float values

    data_pipeline = cluster_data.PrepareClusterData("./data/KSA_human.txt")
    data_pipeline.CancerData = cancer_data
    data_pipeline.replace_with_average()

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


    HierCluster = hierarchical.Hierarchical()
    #add data to hiercluster
    HierCluster.X = X
    HierCluster.labels = labels
    HierCluster.kinaseFeatures = kinaseFeatures
    HierCluster.poorKFeats = poorKFeats
    HierCluster.richKFeats = richKFeats
    HierCluster.labelsPoor = labelsPoor
    HierCluster.labelsRich = labelsRich
    HierCluster.pfile = pfile
    HierCluster.clusterMethod("notpca", 12, pfile)

    #add cluster to hyper geometric class
    comparativeClusterGroups = compare_clusters.CompareCluster(12)
    comparativeClusterGroups.add_cluster(HierCluster.clusters)

    filename = "./data/pickles/randomClusterGroups"
    with open(filename, 'wb+') as f:
        pickle.dump(comparativeClusterGroups.all_clusters, f)

    comparativeClusterGroups.setMainCluster()
    comparativeClusterGroups.create_graph()
    comparativeClusterGroups.get_edge_scores()
    comparativeClusterGroups.display_stats()
    #todo get average of scores
    #!in range 12


    ############Test checking if changing cluster groups #############
    print("cluster groups ", comparativeClusterGroups.all_cluster_nodes[1][1].data)
    
    ###################################################################3


    print("current scores ", scores)
    counter = 0
    for ik in range(len(comparativeClusterGroups.all_cluster_nodes[0])): #12
        for k in range(len(comparativeClusterGroups.all_cluster_nodes[0])): #12
            if indy == 0:
                if ik != k:
                    scores.append(comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[1][k].name])
                else:
                    scores.append(0)

            else:
                scores[counter] += comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[1][k].name]

            counter += 1

    counter = 0
#divide scores by 100 
for i in range(len(scores)):
    scores[i] = scores[i] / 100
   

##pickle data
filename="data/pickles/randomScores"
with open(filename, 'wb+') as f:
    pickle.dump(scores, f)


print("current scores ", scores)


