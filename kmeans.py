from sklearn.cluster import KMeans
import numpy as np
import cluster_data as cd
import jenkspy

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
        x = self.startCluster.get_basal_training_data()
        xCol = []
        #get phosphorylation column data
        for i in range(len(x)):
            xCol.append(x[i][1])

        #reshape data
        xCol = np.array(xCol)
        print(xCol)
        xCol = xCol.reshape((-1,1))
        print(xCol)

        self.kmeans(xCol)

    #jenks algorithm
    def JenksAlgo(self, x):
        breaks = jenkspy.jenks_breaks(x, nb_class=self.clusters)
        print(breaks)

    def runJenks(self):
        x = self.startCluster.get_basal_training_data()

        xCol = []
        #get phosphorylation column data
        for i in range(len(x)):
            xCol.append(x[i][1])


        self.JenksAlgo(xCol)


