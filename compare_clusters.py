import scipy
import pandas as pd
import pickle
import alias
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from scipy.stats import hypergeom
from scipy.special import comb 

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
        self.uniqueKinases = []
        self.all_clusters = []
        self.alias_object = alias.Alias("./data/info_table.csv")
        self.filename = "kinase_class.csv"
        self.family_data = pd.read_csv("./data/kinaseClass.txt", delimiter=",")
        self.family_clusters = {}
        self.all_cluster_nodes = []
        self.overlap = []
        self.k = 0
        self.N = None


    def setMainCluster(self):
        for i in range(len(self.family_data)):
            #!todo check aliases
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



    def add_cluster(self, cluster_group):
        self.all_clusters.append(cluster_group)



    def hyperGeometric(self, overlap, k, N):
        #k total success in M group - overlap
        #n is number of draws(size of 1 correlation cluster)
        #N is total population(all kinases in all families)
        #M is size of one family cluster

        #cumulative hyper geometric distribution
        #x -> current integer in rayman sum
        #M -> Total population of kinases in all clusters (even not in our phospho data)
        #N -> Total population of family cluster i
        #prob = hypergeom.cdf(x, M, n, N)
        M = len(self.family_data)
        #M = len(self.uniqueKinases)
        
        prob = 0
        for i in range(overlap, N):
            #[x] choose [i]
            one = comb(k, i)
            #[M-x] choose [N-i]
            two = comb((M-k), (N-i))
            #[M] choose [N]
            prob += (one * two) / comb(M, N)


        return prob


    def create_graph(self):
        #create nodes for main family tree
        cluster_family_names = list(self.family_clusters.keys())
        cluster_family_vals = list(self.family_clusters.values())
        row = []

        for i in range(len(self.family_clusters)):
            n = Node(cluster_family_names[i])
            n.data = cluster_family_vals[i]
            row.append(n)


        self.all_cluster_nodes.append(row)
        row = []

        #should be two
        #set up unique kinases here
        for i in range(len(self.all_clusters)):
            for j in range(len(self.all_clusters[i])):
                name = "type" + str(i) + "family" + str(j)
                n = Node(name)
                kinaseList = self.all_clusters[i][j]
                n.data = kinaseList
                #check if not in array already
                for kinase in kinaseList:
                    kinase = str(kinase).upper()
                    if kinase not in self.uniqueKinases:
                        self.uniqueKinases.append(kinase)

                row.append(n)

            #add type i group to cluster groups
            self.all_cluster_nodes.append(row)
            row = []

        #pickle kinases for interaction matrix
        filename = "./data/pickles/uniqueKinases"
        if os.path.exists(filename) == False:
            with open(filename, 'wb+') as f:
                pickle.dump(self.uniqueKinases, f)

        #remove kinases not in data
        self.filter_phospho_kinases()

        row = []
        return



    def get_edge_scores(self):
        k = 0
        M = 0
        n = 0
        commonKinases = []
        significantScores = []
        sigNodes = []
        with open("./results/overlapScores.txt", 'w+') as f:
            for i in range(1, len(self.all_cluster_nodes)):
                for kr in range(len(self.all_cluster_nodes)):
                    #computer edge score between all of the ith cluster node and original family cluster node
                    for j in range(len(self.all_cluster_nodes[i])):
                        for s in range(len(self.all_cluster_nodes[kr])):
                            k = 0
                            for kinase in self.all_cluster_nodes[kr][s].data:
                                #if kinase in family is in the cluster group i cluster j add it to overlap
                                if kinase in self.all_cluster_nodes[i][j].data:
                                    commonKinases.append(kinase)
                                    k += 1
                            #before you move to next cluster computer edge score of node
                            #n is length of cluster j
                            n = len(self.all_cluster_nodes[i][j].data)
                            #N is length of family cluster s
                            N = len(self.all_cluster_nodes[kr][s].data)
                            #set edge score of cluster group i and cluster j with family group and cluster s
                            overlap = k
                            score = self.hyperGeometric(overlap, n, N)
                            self.all_cluster_nodes[i][j].edges[self.all_cluster_nodes[kr][s].name] = score
                            if score < .05 and score != 0:
                                significantScores.append(("Edge between " + str(self.all_cluster_nodes[i][j].name) + " and " + str(self.all_cluster_nodes[kr][s].name) + " is " + str(self.all_cluster_nodes[i][j].edges[self.all_cluster_nodes[kr][s].name]) + "\n"))
                                sigNodes.append([self.all_cluster_nodes[i][j].name, self.all_cluster_nodes[kr][s].name])
                            scores = "Edge between " + str(self.all_cluster_nodes[i][j].name) + " and " + str(self.all_cluster_nodes[kr][s].name) + " is " + str(self.all_cluster_nodes[i][j].edges[self.all_cluster_nodes[kr][s].name]) + "\n"
                            f.write(scores)

                        self.overlap.append(commonKinases)
                        commonKinases = []

                        
        #pickle cluster nodes
        with open("./data/pickles/clusterNodes", 'wb+') as fo:
            print("pickling data")
            pickle.dump(self.all_cluster_nodes, fo)


        #significant cluster nodes
        with open("./data/pickles/sigNodes", 'wb+') as f:
            pickle.dump(sigNodes, f)

        with open("./results/significantScores.txt", 'w+') as f:
            for i in range(len(significantScores)):
                f.write(significantScores[i])
    

    #filter out kinases that are not in our phosphorylation data
    def filter_phospho_kinases(self):
        replace_kinases = []
        #in family clusters
        for i in range(len(self.all_cluster_nodes[0])):
            for kinase in self.all_cluster_nodes[0][i].data:
                #uppercase kinase
                kinase = str(kinase).upper()
                if kinase in self.uniqueKinases:
                    replace_kinases.append(kinase)

            self.all_cluster_nodes[0][i].data = replace_kinases
            replace_kinases = []


    def display_stats(self):
        filename = "./results/familyclusters.txt"
        with open(filename, 'w+') as f:
            for i in range(len(self.all_cluster_nodes)):
                for j in range(len(self.all_cluster_nodes[i])):
                    #in each cluster
                    string = "Family " + str(self.all_cluster_nodes[i][j].name) + " Kinase count" +  str(len(self.all_cluster_nodes[i][j].data)) +  " Kinases " +  str(self.all_cluster_nodes[i][j].data) + "\n"
                    f.write(string)


    #check if cluster groups file exists if not create it
    def data_exists_check(self, random=False):
        #look for pickle files
        filename1 = "./data/pickles/clusterGroups"
        if os.path.exists(filename1) and random == False:
            #read in clusters
            with open(filename1, 'rb+') as f:
                self.all_clusters = pickle.load(f)

            return

        else:
            self.get_clusters()


    def get_clusters(self):
        #start hierarchical clustering
        hierCluster = hierarchical.Hierarchical()
        hierCluster.kinaseFile = "./data/KSA_human.txt"
        #hierCluster.dataFile = "./data/BreastCancerData.xlsx"
        hierCluster.clusterMethod("pca", 12)
        hBreastCancerCluster = hierCluster.clusters


        #add hierarchical clustering to all groups
        self.add_cluster(hBreastCancerCluster)


        #get ovarian cancer clustering
        ovHierCluster = hierarchical.Hierarchical()
        ovHierCluster.kinaseFile = "./data/KSA_human.txt"
        ovHierCluster.clusterMethod("pca", 12)
        hOvarianCancerCluster = ovHierCluster.clusters
        self.add_cluster(hOvarianCancerCluster)



        #get colorectal cancer clustering
        '''colHierCluster = hierarchical.Hierarchical()
        colHierCluster.kinaseFile = "./data/KSA_human.txt"
        colHierCluster.clusterMethod("pca", 8)
        colCancerCluster = colHierCluster.clusters
        main.add_cluster(colCancerCluster)'''


        #pickle clusters
        filename = "./data/pickles/clusterGroups"
        with open(filename, 'wb+') as f:
            pickle.dump(self.all_clusters, f)


class Node:
    def __init__(self, name):
        self.data = None
        self.edges = {}
        self.name = name


main = CompareCluster(2)
main.setMainCluster()
main.data_exists_check()
main.create_graph()
main.get_edge_scores()
#main.draw_graph()
main.display_stats()

#hypergf summation cumulative & equal
