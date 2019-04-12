from sklearn.cluster import KMeans

import cluster_data as cd

#parameters 
class Knn_cluster:
    def __init__(self):
        self.clusters = 0
        self.startCluster = cd.ClusterData()

    def knn(self, x):
        knn = KMeans(n_clusters=self.clusters, verbose=1)
        knn.fit(x)
        print(knn.labels_)

    def run_knn(self):
        #grab training data
        x = self.startCluster.get_basal_training_data()
        x = x[:][1]
        self.knn(x)