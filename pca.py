import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from numpy import linalg
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
import cluster_data
import graph
#start pca 

#grab kinase bucket matrixes
def getMatrix(kinase):
    myMatrix = cluster_data.ClusterData(kinase).get_kinase_substrate_matrixes(2)
    return myMatrix


#get SVDs of each kinase bucket and write it to svd txt file
def getSVDdata(kinase):
    kinaseFeature = {}
    matrix = getMatrix(kinase)
    
    with open("svd.txt", 'w+') as f: 
        f.write("Method Singular Value Decomposition(One of the PCA methods)\n")
        f.write("X = U(SIG)V*\n\n")
        f.write("X shape (nxm)\n\n")
        f.write("U shape (nxn)\n\n")
        f.write("(SIG) shape (nx1)\n\n")
        f.write("V shape (nxm)\n\n")

        
        for kinase, data in matrix.items():
            bucket = []
            print("k ", kinase, data)
            for substrate in data.values():
                bucket.append(substrate)

            kinaseFeature[kinase], u, s, vt = getFeatureVector(kinase, bucket)

            f.write("Kinase " + str(kinase) + "\n" + "Singular Vector U \n" + str(u) + "\n" + "Singular Values \n" + str(s) + "\n" + "Singular Vector V (transpose) \n" + str(vt) + "\n\n")
            
    return kinaseFeature

#get principal components of kinase buckets
def getPcaVectors(matrix):
    pca = PCA(n_components=2)
    matrix = StandardScaler().fit_transform(matrix)
    pcs = pca.fit_transform(matrix)
    pcDf = pd.DataFrame(data=pcs, columns=['principal component 1', 'principal component 2'])
    return pcDf


#get kinase feature vector
def getFeatureVector(kinase, matrix):
    #transpose matrix
    matrix = np.transpose(matrix)
    print(matrix.shape)
    #pcs = getPcaVectors(matrix)
    u, s, vt = linalg.svd(matrix, full_matrices=False)
    
    #get first column of V   (6*6) (6*27) ((6*27)->VT)  / (27*27) (27*6) ((27*6) -> VT)  X-> (27*6) V -> (6*27)  VR-> (6 * 1)  X* VR -> (27*6)(6*1) -> (27*1)
    v = np.transpose(vt)
    v = v[:,:1]
    #get scores XVR
    vec = np.matmul(matrix, v)
    return vec, u, s, vt

#list unique kinases and its substrates
def visualizeDataApp(kinase):
    #create kinase and substrate association
    #matrix is dictionary
    matrix = getMatrix(kinase)
    
    #todo only getting two kinases?
    #write to file
    with open('ksa.txt', 'w+') as f:
        for kinase, bucket in matrix.items():
            sub = []

            for substrate, data in bucket.items():
                sub.append(substrate)

            f.write(F'{kinase} {sub}' + "\n")

