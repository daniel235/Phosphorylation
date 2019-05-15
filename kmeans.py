from sklearn.cluster import KMeans
import numpy as np
import cluster_data as cd
from google.cloud import storage
from oauth2client.service_account import ServiceAccountCredentials
import os


#parameters 
class Kmeans_cluster:
    def __init__(self):
        self.clusters = 50
        self.startCluster = cd.ClusterData()

    def kmeans(self, x):
        kmeans = KMeans(n_clusters=self.clusters, verbose=1)
        kmeans.fit(x)
        print(kmeans.cluster_centers_)

    def run_kmeans(self):
        #grab training data
        x = self.startCluster.get_basal_bicor_correlation_matrix()
        shape = np.shape(x)
        #grab indices of upper triangle 
        xIndices = np.triu_indices(shape[0], shape[1])
        
    

    