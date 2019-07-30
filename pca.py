import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from numpy import linalg
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
import cluster_data
import graph
import stats
#start pca 

#grab kinase bucket matrixes
def getMatrix(kinase):
    clusterStructure = cluster_data.ClusterData(kinase)
    myMatrix = clusterStructure.get_kinase_substrate_matrixes(2)
    psiteCount = len(clusterStructure.CancerData[:,1])
    kinaseCount = len(myMatrix.keys())
    tsampleCount = len(clusterStructure.CancerData[1])
    afterStat = stats.Statistics()
    afterStat.set_table(psiteCount, kinaseCount, tsampleCount)
    afterStat.plotTable()
    return myMatrix, clusterStructure.fileName


#get SVDs of each kinase bucket and write it to svd txt file
def getSVDdata(kinase, threshold):
    poorKinaseFeature = {}
    richKinaseFeature = {}
    substrateCount = 0
    kinaseFeature = {}
    matrix, pfile = getMatrix(kinase)
    
    with open("svd.txt", 'w+') as f: 
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


            kinaseFeature[kinase], u, s, vt = getFeatureVector(kinase, bucket, 1)
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
    
    #pcs = getPcaVectors(matrix)
    u, s, vt = linalg.svd(matrix, full_matrices=False)
    
    #get first column of V   (6*6) (6*27) ((6*27)->VT)  / (27*27) (27*6) ((27*6) -> VT)  X-> (27*6) V -> (6*27)  VR-> (6 * 1)  X* VR -> (27*6)(6*1) -> (27*1)
    v = np.transpose(vt)
    #vR
    v = v[:,:dim]
    #get scores XVR
    vec = np.matmul(matrix, v)
    return vec, u, s, vt

