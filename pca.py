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
    afterStat.plotTable()
    #return matrix and name of the file we are using
    return myMatrix, clusterStructure.fileName


#get SVDs of each kinase bucket and write it to svd txt file
def getSVDdata(kinase, threshold, obs=None, pfile=None, matrix=None):
    poorKinaseFeature = {}
    richKinaseFeature = {}
    substrateCount = 0
    kinaseFeature = {}
    
    if obs == None:
        matrix, pfile = getMatrix(kinase)

    
    
    with open("./results/" + str(pfile)[:-5] + "svd.txt", 'w+') as f: 
        f.write("Method Singular Value Decomposition(One of the PCA methods)\n")
        f.write("X = U(SIG)V*\n\n")
        f.write("X shape (nxm)\n\n")
        f.write("U shape (nxn)\n\n")
        f.write("(SIG) shape (nx1)\n\n")
        f.write("V shape (nxm)\n\n")

        
        for kinase, data in matrix.items():
            bucket = []
            substrateCount = len(data.values())
            for substrate in data.values():
                bucket.append(substrate)


            kinaseFeature[kinase], u, s, vt = getFeatureVector(kinase, bucket, 2)
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
    #todo center data
    #todo pickle data
    #?(matrix -> 5000x24)
    #?(u -> 5000x5000)
    #?(s -> 5000x24)
    #?(Vt -> 24x24)
    #?(XxV -> (5000x24)x(24x24))
    #!problem if n is less than m ->  for (nxm matrix)
    #transpose matrix
    matrix = np.transpose(matrix)
    
    #standardize data scale 
    ss = StandardScaler()
    #matrix = ss.fit_transform(matrix)
    #print(matrix)

    #pcs = getPcaVectors(matrix)
    #normalize columns
    mean = 0
    for i in range(len(matrix[0])):
        mean = np.mean(matrix[:,i])
        stdDev = statistics.stdev(matrix[:,i])
        for j in range(len(matrix[:0])):
            matrix[j,i] = (matrix[j,i] - mean) / stdDev

    
    u, s, vt = linalg.svd(matrix, full_matrices=False)
    
    #get variance
    var_explained = np.round(s**2/np.sum(s**2), decimals=3)
    sns.barplot(x=list(range(1,len(var_explained)+1)),
        y=var_explained, color="limegreen")
        
    plt.xlabel("PCS")
    plt.ylabel("Percent Variance Explained")
    plt.savefig('./results/svd_variance.png', dpi=100)
    
    #get first column of V   (6*6) (6*27) ((6*27)->VT)  / (27*27) (27*6) ((27*6) -> VT)  X-> (27*6) V -> (6*27)  VR-> (6 * 1)  X* VR -> (27*6)(6*1) -> (27*1)
    #(6x24) x (24x1)

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