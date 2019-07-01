import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from numpy import linalg
from sklearn.preprocessing import StandardScaler
import pandas as pd
import cluster_data
import graph
#start pca 

#grab kinase matrixes
def getMatrix():
    myMatrix = cluster_data.ClusterData("./data/KSA_human.txt", "./data/BreastCancerData.xlsx").get_kinase_substrate_matrixes(2)
    return myMatrix

#plot histogram count of number of substrates
def getSVDdata():
    matrix = getMatrix()

    with open("svd.txt", 'w+') as f: 
        f.write("Method Singular Value Decomposition(One of the PCA methods)\n")
        f.write("X = U(SIG)V*\n\n")
        f.write("X shape (nxm)\n\n")
        f.write("U shape (nxn)\n\n")
        f.write("(SIG) shape (nx1)\n\n")
        f.write("V shape (nxm)\n\n")

        
        for kinase, data in matrix.items():
            bucket = []
            for substrate in data.values():
                bucket.append(substrate)

            #pc = getPcaVectors(bucket)
            #bucket = np.transpose(np.array(bucket))
            u, s, vt = getSVD(bucket)

            if kinase == "AURKA":
                getVarianceVector(kinase, bucket)

            f.write("Kinase " + str(kinase) + "\n" + "Singular Vector U \n" + str(u) + "\n" + "Singular Values \n" + str(s) + "\n" + "Singular Vector V (transpose) \n" + str(vt) + "\n\n")
            

def getPcaVectors(matrix):
    pca = PCA(n_components=2)
    matrix = StandardScaler().fit_transform(matrix)
    pcs = pca.fit_transform(matrix)
    pcDf = pd.DataFrame(data=pcs, columns=['principal component 1', 'principal component 2'])
    return pcDf
    
def getSVD(matrix):
    u, s, vt = linalg.svd(matrix, full_matrices=False)
    return u, s, vt

def getVarianceVector(kinase, matrix):
    #transpose matrix
    matrix = np.transpose(matrix)
    pcs = getPcaVectors(matrix)
    u, s, vt = getSVD(matrix)
    
    #get first 2 columns of    (6*6) (6*27) ((6*27)->VT)  / (27*27) (27*6) ((27*6) -> VT)  X-> (27*6) V -> (6*27)  VR-> (6 * 1)  X* VR -> (27*6)(6*1) -> (27*1)
    v = np.transpose(vt)
    v = v[:,:1]
    print(v)
    print(matrix)
    vec = np.matmul(matrix, v)
    print(vec)
    for i in range(len(vec)):
        print("vec ", i , vec[i])
        plt.plot(vec[i,0], vec[i, 1], 'o')

    plt.ylabel("pc2")
    plt.xlabel("pc1")
    plt.show()


def compareVectors(pcaMatrix, svdMatrix):
    pass


def visualizeDataApp():
    #create kinase and substrate association
    #matrix is dictionary
    matrix = getMatrix()
    
    #todo only getting two kinases?
    #write to file
    with open('ksa.txt', 'w+') as f:
        for kinase, bucket in matrix.items():
            sub = []

            for substrate, data in bucket.items():
                sub.append(substrate)

            f.write(F'{kinase} {sub}' + "\n")


def printData():
    dataCenter = cluster_data.ClusterData("./data/KSA_human.txt", "./data/BreastCancerData.xlsx")
    matrix = dataCenter.get_kinase_substrate_matrixes(2)
    
    
#printData()
#visualizeDataApp()
getSVDdata()
