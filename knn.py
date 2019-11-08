from sklearn.cluster import KMeans
import numpy as np
import cluster_data as cd
import jenkspy

#parameters 
class Knn_cluster:
    def __init__(self):
        self.clusters = 20
        self.startCluster = cd.ClusterData()

    def knn(self, x):
        knn = KMeans(n_clusters=self.clusters, verbose=1)
        knn.fit(x)
        print(knn.cluster_centers_)

    def run_knn(self):
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

        self.knn(xCol)

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


