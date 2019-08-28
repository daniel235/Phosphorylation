from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial.distance as ssd
from sklearn.cluster import AgglomerativeClustering
from matplotlib import pyplot as plt 
import seaborn as sns
import pickle
import os
import numpy as np

import pca

#todo (correlation clustering  linkage avg)

class Hierarchical:
    def __init__(self):
        self.kinaseFeatures = {}
        self.poorKFeats = {}
        self.richKFeats = {}
        self.kinaseFile = None
        self.dataFile = None
        self.X = []
        self.labels = []
        self.Xpoor = []
        self.labelsPoor = []
        self.Xrich = []
        self.labelsRich = []
        self.pfile = None
        self.k = 0


    #method of clustering/ create distance matrix 
    def clusterMethod(self, method):
        if method == "pca":
            #get kinase svd feature 
            self.kinaseFeatures, self.poorKFeats, self.richKFeats, pfile = pca.getSVDdata(self.kinaseFile, 10)
            for kinase, vector in self.kinaseFeatures.items():
                self.X.append(vector)
                self.labels.append(kinase)

            for kinase, vector in self.poorKFeats.items():
                self.Xpoor.append(vector)
                self.labelsPoor.append(kinase)

            for kinase, vector in self.richKFeats.items():
                self.Xrich.append(vector)
                self.labelsRich.append(kinase)


        self.pfile = pfile
        #plot general kinase clustering
        distMatrix = self.correlationMatrix(self.X, self.labels)

        distArray = ssd.squareform(distMatrix)

        arr = linkage(distArray, method='ward')

        plt.figure()
        dendrogram(arr, labels=self.labels, show_leaf_counts=True, orientation='right', color_threshold=10.0)
        plt.savefig(("./data/results/" + str(pfile)[:-5] + ".jpg"))
        plt.show() 

        #plot poor kinase clustering
        distMatrix = self.correlationMatrix(self.Xpoor, self.labelsPoor)
    
        #condense distance matrix
        distArray = ssd.squareform(distMatrix) 
 
        arr = linkage(distArray, method='ward')
        print(arr)
        #kinase names
        plt.figure()
        dendrogram(arr, labels=self.labelsPoor, show_leaf_counts=True, orientation='right', color_threshold=10.0)
        plt.savefig(("./data/results/" + str(pfile)[:-5] + "poor.jpg"))
        plt.show() 

        #plot rich kinase clustering
        distMatrix = self.correlationMatrix(self.Xrich, self.labelsRich)
    
        #condense distance matrix
        distArray = ssd.squareform(distMatrix) 
 
        arr = linkage(distArray, method='ward')
        #kinase names
        plt.figure()
        dendrogram(arr, labels=self.labelsRich, show_leaf_counts=True, orientation='right', color_threshold=10.0)
        plt.savefig(("./data/results/" + str(self.pfile)[:-5] + "rich.jpg"))
        plt.show() 



    def correlationMatrix(self, matrix, labels):
        #check for pickle
        filename = "./data/pickles/" + str(self.pfile)[:-5] + "correlation"
        if os.path.exists(filename):
            with open(filename, 'rb+') as f:
                correlationMatr = pickle.load(f)
                
        else: 
            "No Correlation Matrix"

        return correlationMatr
        

    def euclidDistance(self, matrix, labels):
        #for i in range(len(a)):
        #get euclid distance of a[i] and b[i]
        distMatrix = []
        row = []
        print("matr", matrix)

        for i in range(len(matrix)):
            a = matrix[i]
            aLabel = labels[i]
            for j in range(len(matrix)):
                b = matrix[j]
                bLabel = labels[j]
                if i == j:
                    row.append(0)
                else:
                    row.append(np.linalg.norm(a-b))

            distMatrix.append(row)
            row = []
    
        return distMatrix


hierCluster = Hierarchical()
hierCluster.kinaseFile = "./data/KSA_human.txt" 
#hierCluster.dataFile = "./data/BreastCancerData.xlsx"
hierCluster.clusterMethod("pca")
