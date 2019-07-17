from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as ssd
from sklearn.cluster import AgglomerativeClustering
from matplotlib import pyplot as plt 
import numpy as np

import pca

class Hierarchical:
    def __init__(self):
        self.kinaseFeatures = {}
        self.kinaseFile = None
        self.dataFile = None
        self.X = []
        self.labels = []
        self.k = 0


    def clusterMethod(self, method):
        if method == "pca":
            #get kinase svd feature 
            self.kinaseFeatures = pca.getSVDdata(self.kinaseFile)
            with open("kFeat.txt", 'w+') as f:
                for kinase, vector in self.kinaseFeatures.items():
                    self.X.append(vector)
                    self.labels.append(kinase)
                    f.write(str(kinase) + "\n\n" + str(vector) + "\n\n")

        distMatrix = self.euclidDistance(self.X)
    
        #condense distance matrix
        distArray = ssd.squareform(distMatrix) 
 
        arr = linkage(distArray, method='single')
        #kinase names
        plt.figure()
        dendrogram(arr, labels=self.labels, show_leaf_counts=True)
        plt.show() 

    def euclidDistance(self, matrix):
        #for i in range(len(a)):
        #get euclid distance of a[i] and b[i]
        distMatrix = []
        row = []
        
        for i in range(len(self.X)):
            a = self.X[i]
            aLabel = self.labels[i]
            for j in range(len(self.X)):
                b = self.X[j]
                bLabel = self.labels[j]
                if i == j:
                    row.append(0)
                else:
                    if(aLabel == 'SYK' and bLabel == 'LCK'):
                        print("syk and lck distance ", np.linalg.norm(a-b))

                    elif aLabel == 'PRKCA' and bLabel == 'MAPK9':
                        print("prkca and mapk9 distance ", np.linalg.norm(a-b))

                    row.append(np.linalg.norm(a-b))

            distMatrix.append(row)
            row = []
    
        return distMatrix
       
        

hierCluster = Hierarchical()
hierCluster.kinaseFile = "./data/KSA_human.txt" 
#hierCluster.dataFile = "./data/BreastCancerData.xlsx"
hierCluster.clusterMethod("pca")
