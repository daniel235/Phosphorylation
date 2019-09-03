import scipy
import pandas as pd
import pickle
import os
import sys

import hierarchical

#compare clusters (for now kmeans)

class CompareCluster:
    '''
        Creates a main cluster that compares to the rest of the cluster groups(hierarchal/kmeans  using correlation)

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
            if self.family_data['Gene'][i] not in self.kinases:
                self.kinases.append(self.family_data['Gene'][i])


            if self.family_data['Classification'][i] not in self.family:
                self.family.append(self.family_data['Classification'][i])

        
        for k in range(len(self.family_data)):
            if self.family_data['Classification'][k] in self.family_clusters:
                self.family_clusters[self.family_data['Classification'][k]].append(self.family_data['Gene'][k])

            else:
                self.family_clusters[self.family_data['Classification'][k]] = [self.family_data['Gene'][k]]

        print(self.family_clusters)
       

    def getNode(self):
        pass


    def setTypes(self):
        pass
        
    def add_cluster(self, cluster_group):
        self.all_clusters.append(cluster_group)


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


#start hierarchical clustering
hierCluster = hierarchical.Hierarchical()
hierCluster.kinaseFile = "./data/KSA_human.txt" 
#hierCluster.dataFile = "./data/BreastCancerData.xlsx"
hierCluster.clusterMethod("pca", 8)
hBreastCancerCluster = hierCluster.clusters


#add hierarchical clustering to all groups
main.add_cluster(hBreastCancerCluster)


#get ovarian cancer clustering
ovHierCluster = hierarchical.Hierarchical()
ovHierCluster.kinaseFile = "./data/KSA_human.txt" 
ovHierCluster.clusterMethod("pca", 8)
hOvarianCancerCluster = ovHierCluster.clusters
main.add_cluster(hOvarianCancerCluster)



#get colorectal cancer clustering
'''colHierCluster = hierarchical.Hierarchical()
colHierCluster.kinaseFile = "./data/KSA_human.txt"
colHierCluster.clusterMethod("pca", 8)
colCancerCluster = colHierCluster.clusters
main.add_cluster(colCancerCluster)'''

print(main.all_clusters)


#pickle clusters
filename = "./data/pickles/clusterGroups"
with open(filename, 'wb+') as f:
    pickle.dump(main.all_clusters, f)
