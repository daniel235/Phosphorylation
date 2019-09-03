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

        Idea - Cytoscape to visualize  sif format

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
        self.all_cluster_nodes = []
        self.k = 0
        self.N = len(self.family_data)


    def setEdge(self, cluster1, cluster2):
        pass


    def setMainCluster(self):
        for i in range(len(self.family_data)):
            #if kinases not in kinase list then add it
            if self.family_data['Gene'][i] not in self.kinases:
                self.kinases.append(self.family_data['Gene'][i])

            #get unique family names
            if self.family_data['Classification'][i] not in self.family:
                self.family.append(self.family_data['Classification'][i])

        #create family dictionary of kinases
        for k in range(len(self.family_data)):
            if self.family_data['Classification'][k] in self.family_clusters:
                self.family_clusters[self.family_data['Classification'][k]].append(self.family_data['Gene'][k])

            else:
                self.family_clusters[self.family_data['Classification'][k]] = [self.family_data['Gene'][k]]

        print(self.family_clusters)
       
        
    def add_cluster(self, cluster_group):
        self.all_clusters.append(cluster_group)
    


    def hyperGeometric(self, k, n, K):
        #observed success x in subgroup
        #total success k in allgroups
        #n is number of draws(size of 1 correlation cluster)
        #N is total population(all kinases in all families)
        success = scipy.misc.comb(k, K)
        failures = scipy.misc.comb(self.N-K,n-k)
        total = scipy.misc.comb(self.N,n)

        prob = (success * failures) / total
        print("probability ", prob)

    def get_edge_scores(self):
        for i in range(1, len(self.all_cluster_nodes)):
            #computer edge score between all of the ith cluster node and original family cluster node
            for j in range(len(self.all_cluster_nodes[i])):
                for k in range(len(self.all_cluster_nodes[0])):
                    self.hyperGeometric()


    def create_graph(self):
        #create nodes for main family tree
        cluster_family_names = list(self.family_clusters.keys())
        cluster_family_vals = list(self.family_clusters.values())
        row = []
        for i in range(len(self.family_clusters)):
            n = Node(cluster_family_names[i])
            n.data = cluster_family_vals[i]
            print("node values ", n.data)
            row.append(n)


        self.all_cluster_nodes.append(row)
        row = []

        #should be two
        for i in range(len(self.all_clusters)):
            for j in range(len(self.all_clusters[i])):
                name = "type" + str(i) + "family" + str(j)
                n = Node(name)
                n.data = self.all_clusters[i][j]
                row.append(n)


        self.all_cluster_nodes.append(row)
        row = []

    def data_exists_check(self):
        #look for pickle files
        filename1 = "./data/pickles/clusterGroups"
        if os.path.exists(filename1):
            return

        else:
            self.get_clusters()

        
    def get_clusters(self):
        #start hierarchical clustering
        hierCluster = hierarchical.Hierarchical()
        hierCluster.kinaseFile = "./data/KSA_human.txt" 
        #hierCluster.dataFile = "./data/BreastCancerData.xlsx"
        hierCluster.clusterMethod("pca", 8)
        hBreastCancerCluster = hierCluster.clusters


        #add hierarchical clustering to all groups
        self.add_cluster(hBreastCancerCluster)


        #get ovarian cancer clustering
        ovHierCluster = hierarchical.Hierarchical()
        ovHierCluster.kinaseFile = "./data/KSA_human.txt" 
        ovHierCluster.clusterMethod("pca", 8)
        hOvarianCancerCluster = ovHierCluster.clusters
        self.add_cluster(hOvarianCancerCluster)



        #get colorectal cancer clustering
        '''colHierCluster = hierarchical.Hierarchical()
        colHierCluster.kinaseFile = "./data/KSA_human.txt"
        colHierCluster.clusterMethod("pca", 8)
        colCancerCluster = colHierCluster.clusters
        main.add_cluster(colCancerCluster)'''

        print(self.all_clusters)


        #pickle clusters
        filename = "./data/pickles/clusterGroups"
        with open(filename, 'wb+') as f:
            pickle.dump(self.all_clusters, f)




class Node:
    def __init__(self, name):
        self.data = None
        self.edges = {}
        self.name = None




main = CompareCluster(2)
main.setMainCluster()
main.data_exists_check()
main.create_graph()