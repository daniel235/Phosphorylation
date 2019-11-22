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
        self.clusters = []


    #method of clustering/ create distance matrix 
    def clusterMethod(self, method, pfile=None):
        if method == "pca":
            #get kinase svd feature 
            self.kinaseFeatures, self.poorKFeats, self.richKFeats, pfile = pca.getSVDdata("./data/KSA_human.txt", 10)
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



    def start_hierarchical_clustering(self, cluster_len, show_plots=False):
        #plot general kinase clustering
        distMatrix = self.correlationMatrix(self.kinaseFeatures, 'all')
        
        distArray = ssd.squareform(distMatrix, checks=False)

        arr = linkage(distArray, method='ward')

        #agglomerative
        cluster = AgglomerativeClustering(n_clusters=cluster_len, affinity='precomputed', linkage='complete')
        print(cluster.fit_predict(distMatrix))

        #initialize empty arrays for each numbered cluster ex: 1, 2, 3
        for i in range(cluster_len):
            self.clusters.append([])

        #create cluster group
        keys = list(self.kinaseFeatures.keys())
        for i in range(len(cluster.labels_)):
            #match to kinase label
            #if kinase not in cluster already add it
            if keys[i] not in self.clusters[cluster.labels_[i]]:
                self.clusters[cluster.labels_[i]].append(keys[i])

        

        plt.figure()
        plot_file = str(self.pfile).replace(".txtls", '')
        dendrogram(arr, labels=self.labels, show_leaf_counts=True, orientation='right', color_threshold=10.0)
        plt.savefig(("./data/results/" + plot_file + ".png"))
        if show_plots:
            plt.show() 

        plt.close()
        #plot poor kinase clustering
        distMatrix = self.correlationMatrix(self.poorKFeats, 'poor')
    
        #condense distance matrix
        distArray = ssd.squareform(distMatrix) 
 
        arr = linkage(distArray, method='ward')
        #kinase names
        plt.figure()
        dendrogram(arr, labels=self.labelsPoor, show_leaf_counts=True, orientation='right', color_threshold=10.0)
        plt.savefig(("./data/results/" + plot_file + "poor.png"))
        if show_plots:
            plt.show() 
        plt.close()
        #plot rich kinase clustering
        distMatrix = self.correlationMatrix(self.richKFeats, 'rich')

        #condense distance matrix
        distArray = ssd.squareform(distMatrix) 
 
        arr = linkage(distArray, method='ward')
        #kinase names
        plt.figure()
        dendrogram(arr, labels=self.labelsRich, show_leaf_counts=True, orientation='right', color_threshold=10.0)
        plt.savefig(("./data/results/" + plot_file + "rich.png"))
        if show_plots:
            plt.show() 
        plt.close()


    def correlationMatrix(self, kfeats, type):
        correlationMatr = None
        #check for pickle
        filename = "./data/pickles/" + str(self.pfile)[:-5] + str(type) + "correlation"
        if os.path.exists(filename):
            with open(filename, 'rb+') as f:
                correlationMatr = pickle.load(f)
                
        else: 
            correlationMatr = self.correlation(kfeats, filename)

        #get euclid distance from correlation matrix
        euclidMatr = []
        row = []
        for i in range(len(correlationMatr)):
            for j in range(len(correlationMatr)):
                if i != j:
                    row.append(np.linalg.norm(correlationMatr[i]-correlationMatr[j]))

                else:
                    row.append(0)


            euclidMatr.append(row)
            row = []


        return euclidMatr

    #?getting correlation from projected data
    def correlation(self, kfeatures, filename):
        fileName = filename
        #fileName = './data/pickles/' + str(self.pfile)[:-5] + 'correlation'
        if os.path.exists(fileName):
            f = open(fileName, 'rb+')
            corrmatr = pickle.load(f)
            f.close()
            return corrmatr


        corrmatr = np.zeros((len(kfeatures.keys()), len(kfeatures.keys())))
        a = 0
        b = 0
        temp1 = []
        temp2 = []
        vec1 = []
        vec2 = []
        second_vec1 = []
        second_vec2 = []
        for kinase, vector in kfeatures.items():
            b = 0
            for coord in vector:
                vec1.append(coord[0])
                vec2.append(coord[1])


            for kinase2, vector2 in kfeatures.items():
                #combine vectors
                for coord in vector2:
                    second_vec1.append(coord[0])
                    second_vec2.append(coord[1])
                    
                #combine matrices
                combined_length = len(vec1) + len(second_vec1)
    
                for i in range(combined_length):
                    if i >= len(vec1):
                        temp1.append(second_vec1[i-len(vec1)])
                        temp2.append(second_vec2[i-len(vec1)])
                    else:
                        temp1.append(vec1[i])
                        temp2.append(vec2[i])
                
                if a != b:
                    corrmatr[a, b] = np.corrcoef(temp1, temp2)[0][1]
                else:
                    corrmatr[a,b] = 0

                second_vec1 = []
                second_vec2 = []
                temp1 = []
                temp2 = []
                b += 1

            a += 1
            vec1 = []
            vec2 = []
            temp1 = []
            temp2 = []
      
        #pickle corrmatr
        fileName = './data/pickles/' + str(self.pfile)[:-5] + 'correlation'
        f = open(fileName, 'wb+')
        pickle.dump(corrmatr, f)
        f.close()
        return corrmatr
        

    def euclidDistance(self, matrix, labels):
        #for i in range(len(a)):
        #get euclid distance of a[i] and b[i]
        distMatrix = []
        row = []

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


    def get_dendrogram_clusters(self, linkage):
        pass

'''
hierCluster = Hierarchical()
hierCluster.kinaseFile = "./data/KSA_human.txt" 
#hierCluster.dataFile = "./data/BreastCancerData.xlsx"
hierCluster.clusterMethod("pca", 12)
print("done")'''
