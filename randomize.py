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
im = np.zeros(shape=(12,12))

for indy in range(50):
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
    HierCluster.clusterMethod("notpca", pfile)
    

    #add cluster to hyper geometric class
    comparativeClusterGroups = compare_clusters.CompareCluster()
    comparativeClusterGroups.add_cluster(HierCluster.clusters)

    filename = "./data/pickles/randomClusterGroups"
    with open(filename, 'wb+') as f:
        pickle.dump(comparativeClusterGroups.all_clusters, f)

    comparativeClusterGroups.setMainCluster()
    comparativeClusterGroups.create_graph()
    comparativeClusterGroups.get_edge_scores(random=True)
    comparativeClusterGroups.display_stats()
    #todo get average of scores
    #!in range 12

    #?start random interaction
    matrix = interactionMatrix.InteractionMatrix(comparativeClusterGroups.all_clusters)
    im = np.add(im, matrix.run_interaction())
    

    counter = 0
    for ik in range(len(comparativeClusterGroups.all_cluster_nodes[0])): #12
        for k in range(len(comparativeClusterGroups.all_cluster_nodes[0])): #12
            if indy == 0:
                if ik != k:
                    scores.append(comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[0][k].name])
                else:
                    scores.append(0)

            else:
                scores[counter] += comparativeClusterGroups.all_cluster_nodes[1][ik].edges[comparativeClusterGroups.all_cluster_nodes[0][k].name]

            counter += 1

print("finished run")
    
#divide scores by 100 
for i in range(len(scores)):
    scores[i] = scores[i] / 50


#averaging interaction matrix
im = np.floor_divide(im, 50)
matrix.save_matrix(im, matrix.family)



##pickle data
filename="data/pickles/randomScores"
with open(filename, 'wb+') as f:
    pickle.dump(scores, f)


print("current scores ", scores)

def get_sig_scores(obj, cg):
    obj = np.array(obj)
    

    ##remove zeros from low list
    for k in range(len(obj)):
        if obj[k] == 0:
            #add 5 to it
            obj[k] = 5
    

    o = obj.argsort()
    indexes = []
    lowScores = []
    for p in range(5):
        print(obj[o[p]])
        lowScores.append(obj[o[p]])
        indexes.append(o[p])
        

    print(indexes)
    sigNodes = []
    for j in range(len(indexes)):
        index1 = int(math.floor(indexes[j] / 13))
        index2 = int(math.floor(indexes[j] % 13)) 
        #!todo not going to find it
        sigNodes.append([cg[1][index1].name, cg[0][index2].name])

    #write sig nodes to pickle
    with open("./data/pickles/randomSignificantNodes", 'wb+') as f:
        pickle.dump(sigNodes, f)


    print(sigNodes)

with open("./data/pickles/randomScores", 'rb+') as f:
    obj = pickle.load(f)

with open("./data/pickles/clusterNodes", 'rb+') as f:
    cgNodes = pickle.load(f)


def distro_scores(scores):
    sns.distplot(scores)


get_sig_scores(obj, cgNodes)
#distro_scores(obj)

