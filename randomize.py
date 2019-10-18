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


scores = []


cancer_data = np.array(pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object))
print(cancer_data)

#clean data
cleanData = clean.cleanMatrix()
cleanData.data = cancer_data
cleanData.set_gene_site_column([0,1], True)
cleanData.omit_columns([16,17,18])
cleanData.clean_rows()

cancer_data = cleanData.data


for i in range(100):
    print("run ", i)
    #insert random data into cells
    print(cancer_data[0])
    for i in range(len(cancer_data)):
        for j in range(1, len(cancer_data[i])):
            cancer_data[i][j] = random.uniform(-2, 3)
            #todo need to get float values
            
    print(cancer_data[0])

    data_pipeline = cluster_data.PrepareClusterData("./data/KSA_human.txt")
    print("got past pipeline")
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

    print("before hier")
    HierCluster = hierarchical.Hierarchical()
    print("after hier")
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
    for ik in range(len(comparativeClusterGroups.all_cluster_nodes[0])):
        for k in range(len(comparativeClusterGroups.all_cluster_nodes[0])):
            if len(scores) < len(comparativeClusterGroups.uniqueKinases):
                scores.append(comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[0][k]]])

            else:
                scores[k] += comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[ik][]]


    print("current scores ", scores)


