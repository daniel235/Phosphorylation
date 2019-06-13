import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import cluster_data
import graph
#start pca 

#grab kinase matrixes
def getMatrix():
    myMatrix, subLength = cluster_data.ClusterData().get_kinase_substrate_matrixes(1)
    return myMatrix, subLength

#plot histogram count of number of substrates
def plotHistogram(n):
    matrix, subLength = getMatrix()
    weak, rich = graph.Graph().createKinaseClassHistograms(n)
    bucket = []
    #plt.hist(rich, bins='auto')
    #plt.show()
    
    principalComponents = []
    for kinase in matrix.values():
        bucket = []
        for substrate in kinase.values():
            print(len(substrate))
            bucket.append(substrate)

        
        pc = getPcaVectors(bucket)
        principalComponents.append(pc)
        bucket = []

    #print(principalComponents, "\n")

def getPcaVectors(matrix):
    pca = PCA(n_components=1)
    matrix = StandardScaler().fit_transform(matrix)
    pcs = pca.fit_transform(matrix)
    pcDf = pd.DataFrame(data=pcs, columns=['principal component 1'])
    return pcDf
    
def visualizeDataApp():
    #create kinase and substrate association
    #matrix is dictionary
    matrix, s = getMatrix()

    #todo only getting two kinases?
    #write to file
    with open('ksa.txt', 'w+') as f:
        for kinase, bucket in matrix.items():
            sub = []
            for substrate in bucket:
                sub.append(substrate)
            
            f.write(F'{kinase} {sub}')


def printData():
    dataframe = cluster_data.ClusterData().breastCancerData


printData()
