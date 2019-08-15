from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import cluster_data as cd
import pca
import os
import pickle

#kmediods

#parameters 
class Kmeans_cluster:
    def __init__(self, kinaseFile):
        self.cluster_name = ""
        self.nclusters = 10
        self.X = []
        self.labels = []
        self.kinaseFeatures = {}
        self.poorKinaseFeatures = {}
        self.richKinaseFeatures = {}
        self.clusters = {}
        self.correlationMatrix = None
        self.kinaseFile = kinaseFile


    def run_kmeans(self):
        kmeans = KMeans(n_clusters=self.nclusters, verbose=1)
        kmeans.fit(self.X)
        kmeans.labels_ = self.labels
        print(kmeans.predict(self.X))
        one = []
        kpredictions = kmeans.predict(self.X)
        count = 0
        for k in kpredictions: 
            if k not in self.clusters:
                self.clusters[k] = {}
            
            if self.labels[count] not in self.clusters[k]:
                self.clusters[k][self.labels[count]] = 1

            else:
                self.clusters[k][self.labels[count]] += 1
            
            count += 1

        filename = "./results/" + str(self.cluster_name)[:-5] + "kmeanscluster.txt"
        with open(filename, 'w+') as f:
            for k, kinase in sorted(self.clusters.items()):
                sentence = "\nCluster " + str(k) + "\n"
                f.write(sentence)
                for kinase, data in kinase.items():
                    sentence = str(kinase) + " " + str(data) + " "
                    f.write(sentence)


    def prepareKmeansCluster(self, method):
        if method == "pca":
            #check for data
            filename = input("Enter File Name \n")
            files = "./data/pickles/kinaseFeatures" + str(filename)[:-5]
            
            if os.path.exists(files) == False:
                f = open(files, 'wb+')
                self.kinaseFeatures, self.poorKinaseFeatures, self.richKinaseFeatures, pfile = pca.getSVDdata(self.kinaseFile, 10)
                #pickle data
                pickle_dict = {'kf': self.kinaseFeatures, 'pkf': self.poorKinaseFeatures, 'rkf': self.richKinaseFeatures}
                pickle.dump(pickle_dict, file=f)
                f.close()
            else:
                f = open(files, 'rb+')
                kinase_feat_dict = pickle.load(f, encoding='latin1')
                self.kinaseFeatures = kinase_feat_dict['kf']
                self.poorKinaseFeatures = kinase_feat_dict['pkf']
                self.richKinaseFeatures = kinase_feat_dict['rkf']

            
            if filename not in os.listdir("./data"):
                raise ValueError('No File by this name')
               
            else:
                self.cluster_name = filename
                
            corrmatr = self.correlation(self.kinaseFeatures)
            self.correlationMatrix = corrmatr
            for kinase, vector in self.kinaseFeatures.items():
                print("vec ", vector)
                for x in vector:
                    self.X.append(x)
                    self.labels.append(kinase)
                    

        self.X = np.array(self.X)
        self.kmeansPlot()
        return


    def kmeansPlot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        print("before log ", self.X[1:5])
        self.logValues()
        print("after log " , self.X[1:5])

        #for i in range(len(self.X)):
        fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        '''xs = self.X[i][:,0]
        ys = self.X[i][:,1]
        zs = self.X[i][:,2]
        xs2 = self.X[j][:,0]
        ys2 = self.X[j][:,1]
        zs2 = self.X[j][:,2]
        ax.scatter(xs, ys, zs, 'o', label=self.labels[i])
        ax.scatter(xs2, ys2, zs2, 'o', label=self.labels[j])
        '''
        xbucket = []
        ybucket = []
        currentlabel = self.labels[0]
        for i in range(len(self.labels)):
            if self.labels[i] != currentlabel:
                plt.plot(xbucket, ybucket, 'o', label=currentlabel)
                xbucket = []
                ybucket = []
                currentlabel = self.labels[i]

            else:
                xbucket.append(self.X[i,0])
                ybucket.append(self.X[i,1])
        
            
        plt.legend()
        filename = "./results/" + str(self.cluster_name)[:-5] + "kmeans.png"
        plt.savefig(filename)
        plt.show()
        
    
    def correlation(self, kfeatures):
        fileName = './data/pickles/' + str(self.cluster_name)[:-5] + 'correlation'
        if os.path.exists(fileName):
            f = open(fileName, 'rb+')
            corrmatr = pickle.load(f)
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
                    print("temp1 ", temp1)
                    corrmatr[a, b] = np.corrcoef(temp1, temp2)[0][1]
                else:
                    corrmatr[a,b] = 1

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
        
        print("corr matr 2, 2 ", corrmatr[2,2])
        #pickle corrmatr
        fileName = './data/pickles/' + str(self.cluster_name)[:-5] + 'correlation'
        f = open(fileName, 'wb+')
        pickle.dump(corrmatr, f)
        return corrmatr


    def logValues(self):
        for i in range(len(self.X)):
            if self.X[i][0] <= 0:   
                self.X[i][0] = -1 * math.log2(abs(self.X[i][0]))
                self.X[i][1] = -1 * math.log2(abs(self.X[i][1]))

            else:
                self.X[i][0] = math.log2(abs(self.X[i][0]))
                self.X[i][1] = math.log2(abs(self.X[i][1]))
            

    def heatmapCluster(self):
        cg = sns.clustermap(self.correlationMatrix, cmap="YlGnBu", linewidths=0.1)
        plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        plt.show()
        

Kmeans = Kmeans_cluster("./data/KSA_human.txt")
Kmeans.prepareKmeansCluster("pca")
Kmeans.run_kmeans()
Kmeans.heatmapCluster()
        
    

    