from sklearn.cluster import KMeans
import numpy as np
import cluster_data as cd
import pca
import os


#kmediods

#parameters 
class Kmeans_cluster:
    def __init__(self, kinaseFile):
        self.clusters = 50
        self.startCluster = cd.ClusterData(kinaseFile)
        self.X = []
        self.labels = []
        self.kinaseFeatures = {}
        self.poorKinaseFeatures = {}
        self.richKinaseFeatures = {}
        self.kinaseFile = kinaseFile

    def kmeans(self, x):
        kmeans = KMeans(n_clusters=self.clusters, verbose=1)
        kmeans.fit(x)
        print(kmeans.cluster_centers_)


    def run_kmeans(self):
        #grab training data
        x = None
        shape = np.shape(x)
        #grab indices of upper triangle 
        xIndices = np.triu_indices(shape[0], shape[1])

    
    def graph_k(self):
        pass

    def kmeansCluster(self, method):
        if method == "pca":
            self.kinaseFeatures, pfile = pca.getSVDdata(self.kinaseFile)
            for kinase, vector in self.kinaseFeatures.items():
                self.X.append(vector)
                self.labels.append(kinase)

        print(self.X)
        return

Kmeans = Kmeans_cluster("./data/KSA_human.txt")
Kmeans.kmeansCluster("pca")

        
    

    