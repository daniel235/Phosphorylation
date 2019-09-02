import scipy
import pandas as pd
import os
import sys

import kmeans

#compare clusters (for now kmeans)

class CompareCluster:
    '''
        Creates a main cluster that compares to the rest of the clusters

        O O O O O O     nodes(nclusters)
        |/|\|\| | | |     edges (# of matching kinases)  #compare one cluster with all clusters
        O O O O O O     nodes(nclusters)

        Arguments: 
            Set of arrays are clusters ['AURKA', 'LSK', 'PC1']
            ['UMB', 'CDK1', 'CDK2', 'CDK5', 'HEM1']
            ['SYK', 'LCK']

        Idea - Cytoscape to visualize

    '''
    def __init__(self, nclusters):
        self.accuracy = None
        self.types = []
        self.nclusters = nclusters
        self.kinases = []
        self.family = []
        self.all_clusters = []
        self.filename = "kinase_class.csv"
        self.family_data = pd.read_csv("./data/kinaseClass.txt", delimiter=",")
        self.family_clusters = {}

    def setEdge(self, cluster1, cluster2):
        pass


    def setMainCluster(self):
        for i in range(len(self.family_data)):
            #if kinases not in kinase list then add it
            if self.family_data[i]['Gene'] not in self.kinases:
                self.kinases.append(self.family_data[i]['Gene'])


            if self.family_data[i]['Classification'] not in self.family:
                self.family.append(self.family_data[i]['Classification'])

        
        for k in range(len(self.family_data)):
            if self.family_data[k]['Classification'] in self.family_clusters:
                self.family_clusters[self.family_data[k]['Classification']].append(self.family_data[k]['Gene'])

            else:
                self.family_clusters[self.family_data[k][1]] = [self.family_data[k][0]]

        print(self.family_clusters)

    def getNode(self):
        pass


    def setTypes(self):
        pass
        
    def add_cluster(self, cluster):
        self.all_clusters.append(cluster)


    def hyperGeometric(self, k, N, n, K):
        #observed success x in sub
        #total success k in sub
        #n is number of draws
        #N is total population
        success = scipy.misc.comb(k, K)
        failures = scipy.misc.comb(N-K,n-k)
        total = scipy.misc.comb(N,n)

        prob = (success * failures) / total
        print("probability ", prob)



main = CompareCluster(2)
main.setMainCluster()


#add cluster to object


#pd.read_csv(r"./kinaseClass.txt")
