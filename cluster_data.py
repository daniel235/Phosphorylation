import numpy as np
from sklearn.utils import shuffle
import pipe_line as pipe
import math
import pandas as pd
from scipy.stats import stats
from astropy import stats
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
            fileObject = open(fileName, 'r')
            pearson = pickle.load(fileObject)
            fileObject.close()
        else:
            #pearson correlation
            pearson = np.triu(np.corrcoef(basal_data))

            #write bcors to pickle
            fileObject = open(fileName, 'wb')

            pickle.dump(pearson, fileObject)

            #close file
            fileObject.close()
        
        
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



        


        


