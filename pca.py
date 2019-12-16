import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from numpy import linalg
from scipy import linalg
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
import pickle
import pandas as pd
import seaborn as sns
import statistics
import os
import cluster_data
import graph
import stats
import math
#start pca 

Kinase_variance_vectors = {}
pfile = None

#grab kinase bucket matrixes
def getMatrix(kinase, test=False):
    #clean data first
    #!check this code
    clusterStructure = cluster_data.PrepareClusterData(kinase)
    clusterStructure.clean_data()
    #grab substrate matrixes with minimum requried substrates n
    myMatrix = clusterStructure.get_kinase_substrate_matrixes(2)
    #how many psites are in data after it is cleaned
    psiteCount = len(clusterStructure.CancerData[:,1])
    #how many kinases are matched to the phosphorylation data
    kinaseCount = len(myMatrix.keys())
    #how many tumor samples are in the data after it is cleaned
    tsampleCount = len(clusterStructure.CancerData[1]) - 1
    #get statistics (how many psites/kinases are in this filtered data)
    afterStat = stats.Statistics() 
    afterStat.set_table(psiteCount, kinaseCount, tsampleCount)
    #afterStat.plotTable()
    #return matrix and name of the file we are using
    return myMatrix, clusterStructure.fileName


#get SVDs of each kinase bucket and write it to svd txt file
def getSVDdata(kinase, threshold, obs=None, pfile=None, matrix=None):
    #*necessary data structures
    poorKinaseFeature = {}
    richKinaseFeature = {}
    substrateCount = 0
    kinaseFeature = {}
    
    if obs == None:
        matrix, pfile = getMatrix(kinase)

    #write principal components to svd text file
    with open("./results/" + str(pfile)[:-5] + "svd.txt", 'w+') as f: 
        f.write("Method Singular Value Decomposition(One of the PCA methods)\n")
        
        for kinase, data in matrix.items():
            bucket = []
            substrateCount = len(data.values())
            for substrate in data.values():
                bucket.append(substrate)


            kinaseFeature[kinase], u, s, vt = getFeatureVector(kinase, bucket, 2)
            #separate out poor and rich kinases
            if substrateCount > threshold:
                richKinaseFeature[kinase] = kinaseFeature[kinase]
            else:
                poorKinaseFeature[kinase] = kinaseFeature[kinase]

            f.write("Kinase " + str(kinase) + "\n" + "Singular Vector U \n" + str(u) + "\n" + "Singular Values \n" + str(s) + "\n" + "Singular Vector V (transpose) \n" + str(vt) + "\n\n")
            
    
    return kinaseFeature, poorKinaseFeature, richKinaseFeature, pfile

#get principal components of kinase buckets
def getPcaVectors(matrix):
    pca = PCA(n_components=2)
    matrix = StandardScaler().fit_transform(matrix)
    pcs = pca.fit_transform(matrix)
    pcDf = pd.DataFrame(data=pcs, columns=['principal component 1', 'principal component 2'])
    return pcDf


#get kinase feature vector
def getFeatureVector(kinase, matrix, dim):
    #transpose matrix
    matrix = np.transpose(matrix)

    #normalize columns
    for i in range(len(matrix[0])):
        mean = np.mean(matrix[:,i])
        stdDev = statistics.stdev(matrix[:,i])
        for j in range(len(matrix[:0])):
            print("in get feat vec ", j)
            matrix[j,i] = (matrix[j,i] - mean) / stdDev

    
    u, s, vt = linalg.svd(matrix, full_matrices=False)
    
    #get variance
    var_explained = np.round(s**2/np.sum(s**2), decimals=3)
    sns.barplot(x=list(range(1,len(var_explained)+1)),
        y=var_explained, color="limegreen")
        
    plt.xlabel("PCS")
    plt.ylabel("Percent Variance Explained")
    plt.savefig('./results/svd_variance.png', dpi=100)
    
    #transpose matrix 
    v = np.transpose(vt)
    #vR
    v = v[:,:dim]
    #get scores XVR
    vec = np.matmul(matrix, v)
    

    if os.path.exists('./data/pickles/xvr') == False:
        dbfile = open('./data/pickles/xvr', 'ab')

        pickle.dump(vec, dbfile)
        dbfile.close()

    variance_vector(kinase, matrix)
        
    return vec, u, s, vt


def variance_vector(kinase, matrix):
    #transpose matrix back to original form
    matrix = np.transpose(matrix)
    u, s, vt = linalg.svd(matrix, full_matrices=True)
    Kinase_variance_vectors[kinase] = vt[0]
    #pickle variance vector
    filename = "./data/pickles/" + str(pfile)[:3] + "variance_vector"
    with open(filename, 'wb+') as f:
        pickle.dump(Kinase_variance_vectors, f)