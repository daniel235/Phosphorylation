import numpy as np
from sklearn.utils import shuffle
import pipe_line as pipe
import math
import pandas as pd
from scipy.stats import stats
from astropy import stats
from sklearn.decomposition import PCA
import pickle
import os


class ClusterData:
    '''This class creates the prepared data for the knn and also hierarchy clustering
    functions.  

    All functions are separated by underscores

    returns datasets for basal and luminal types.

    Also creates the y data to test against data.

    Automatically separates data into train and test set

    
    
    Returns:
        prepared data for clustering algorithms as class properties -> (self.trainBasal)
    '''

    def __init__(self):
        self.basal_data = None
        self.pipe_object = pipe.Pipe_line()
        self.basal_pVector = []
        self.luminal_pVector = []
        self.trainBasal, self.trainLuminal = self.get_training_data()
        self.testBasal = None
        self.testLuminal = None
        self.strongKinase = []
        self.weakKinase = []
        self.kinaseData = np.array(pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t"))
        self.phosphorylationData = np.array(pd.read_csv("./data/phosphorylation_data.txt", delimiter="\t"))

    #grab training data from phosphosite database
    def get_training_data(self):
        phosphosite_data = np.array(self.pipe_object.phosphorylation)
        basal_dataPoints = []
        luminal_dataPoints = []
        
        #grab dimension of data
        for count in range(len(phosphosite_data)):
            p = phosphosite_data[count]
            temp = []
            temp2 = []
            for i in range(1, 19):
                if i <= 9:
                    #adding name of psite and phosphorylation data
                    basal_dataPoints.append([p[0],p[i]])
                    temp.append(p[i])
                else:
                    luminal_dataPoints.append([p[0], p[i]])
                    temp2.append(p[i])
 
            self.basal_pVector.append(temp)
            self.luminal_pVector.append(temp2)

        #shuffle data 
        basal_dataPoints, luminal_dataPoints = shuffle(basal_dataPoints, luminal_dataPoints)
        trainTestRatio = math.floor(len(basal_dataPoints) * .7)
        self.trainBasal = basal_dataPoints[0:trainTestRatio]
        self.testBasal = basal_dataPoints[trainTestRatio:]
        self.trainLuminal = luminal_dataPoints[0:trainTestRatio]
        self.testLuminal = luminal_dataPoints[trainTestRatio:]

        return basal_dataPoints, luminal_dataPoints

    ##########  Basal functions  ##############
    def get_basal_training_data(self):
        return self.trainBasal

    def get_basal_testing_data(self):
        return self.testBasal

    def get_y_basal_data(self):
        #data with a y
        y_basal_data = []
        #find ksubs
        for psite in self.trainBasal:
            index = self.pipe_object.find_kinase(psite)

    #correlation kmeans


    ########### Luminal Functions #############
    #leave out test set
    def get_luminal_training_data(self):
        return self.trainLuminal


    def get_luminal_testing_data(self):
        pass

    def get_y_luminal_data(self):
        pass
    

    #!bicor function data
    def get_basal_bicor_correlation_matrix(self):
        basal_data = self.basal_pVector
        #check for bcor value
        fileName = "./pickles/bcorPickle"
        if(os.stat(fileName).st_size != 0):
            with open(fileName, 'rb') as f:
                pearson = pickle.load(f)
           
        else:
            #pearson correlation
            pearson = np.triu(np.corrcoef(basal_data))

            #write bcors to pickle
            fileObject = open(fileName, 'wb')

            pickle.dump(pearson, fileObject)

            #close file
            fileObject.cl23ose()
        
        
        print(np.shape(pearson))
        return pearson

    #!bicor function data
    def get_luminal_correlation_matrix(self):
        luminal_data = self.luminal_pVector

        pearson = np.triu(np.corrcoef(luminal_data))

        return pearson

    #set number of substrates to put kinases in rich/poor class
    def set_arbitrary_kinase_class(self, n):
        kinaseDict = {}
        names = []
        kinaseData = pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t")
        
        #set numbers of 
        for i in range(len(kinaseData)):
            keys = kinaseData["Kinase"][i]
            if keys in kinaseDict: 
                kinaseDict[keys] += 1

            else:
                kinaseDict[keys] = 1
                names.append(keys)

    
        print(kinaseDict)

        #separate kinases
        for i in range(len(kinaseDict)):
            if kinaseDict[names[i]] > n:
                self.strongKinase.append(names[i])
                print("strong ", names[i])
            else:
                self.weakKinase.append(names[i])
                print("weak ", names[i])


    def pca(self, corrMatr):
        #svd(single value decomposition) of correlation matrix
        u, s, vh = np.linalg.svd(corrMatr)
        print(u.shape)
        print(s.shape)
        print(vh.shape)
        print(s)

        #get rank cutoff for data dimension
        r = 2

        pca = PCA(n_components=r)

        pcs = pca.fit_transform(corrMatr)
        pcDf = pd.DataFrame(data=pcs, columns=['principal component 1', 'principal component 2'])
        print(pcDf)

        return u, s, vh

    def count_substrates(self, kinase, ordered=False):
        start = False
        subCount = 0
        indexRange = []
        for i in range(len(self.kinaseData)):
            if self.kinaseData[i][0] == kinase:
                start = True

            #this should stop loop once kinase name is passed (since its alphabetical)
            elif(start and ordered and self.kinaseData[i][0] != kinase):
                break
                
    #?grab substrates of kinase passed in
    def grab_substrates(self, kinase, kinaseFileOrdered=False, PhosDataOrdered=False):
        substrate_matrix = []
        substrate_names = []
        start = False
    
        for i in range(len(self.kinaseData)):
            if self.kinaseData[i][0] == kinase:
                #?logic controllers
                mid = int(len(self.phosphorylationData) / 2)
                current = mid
                start = True
                prevMid = 0
                substrate = self.kinaseData[i][1]
                print("sub ", substrate)
            
                #find substrate in phosphorylation data
                #binary search the data
                if PhosDataOrdered:
                    mins = 0
                    while(mid > 0 and prevMid != mid):
                        #if not found
                        if prevMid == mid:
                            break

                        #substrate at current index
                        index = self.phosphorylationData[int(current)][0]
                        
                        #didn't find substrate
                        if substrate != index:
                            prevMid = mid
                            mid = int(mid / 2)
                            
                            #move to top portion of list
                            if substrate > index:
                                #set index to half of mid + mid
                                mins = current
                                current = int(mid + mins)

                            #move to bottom portion of list
                            else:
                                current = int(mid + mins)
                
                        #found substrate
                        else:
                            #add substrate to names and add it's data row to matrix
                            substrate_names.append(substrate)
                            data = []
                            for i in range(1, len(self.phosphorylationData[current])):
                                data.append(self.phosphorylationData[current][i])

                            substrate_matrix.append(data)
                            break
                
                #brute force search the data
                else:
                    pass


                
            #this should stop loop once kinase name is passed (since its alphabetical)
            elif(start and kinaseFileOrdered and self.kinaseData[i][0] != kinase):
                break

        return substrate_names, substrate_matrix


                

    def get_kinase_substrate_matrixes(self, threshold):
        kinase_matrixes = []
        substrate_threshold = threshold
        kinases = list(set(self.kinaseData[:,0]))
        for i in range(len(kinases)):
            self.grab_substrates(kinases[i], True, True)

        print(kinases)



        


