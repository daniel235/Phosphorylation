import numpy as np
from sklearn.utils import shuffle
import pipe_line as pipe


class ClusterData:
    def __init__(self):
        self.basal_data = None
        self.pipe_object = pipe.Pipe_line()
        self.trainBasal, self.trainLuminal = self.get_training_data()
        self.testBasal = None
        self.testLuminal = None
        
    #grab training data from phosphosite database
    def get_training_data(self):
        phosphosite_data = np.array(self.pipe_object.phosphorylation)
        basal_dataPoints = []
        luminal_dataPoints = []
        
        #grab dimension of data
        for count in range(len(phosphosite_data)):
            p = phosphosite_data[count]
            for i in range(1, 19):
                if i <= 9:
                    #adding name of psite and phosphorylation data
                    basal_dataPoints.append([p[0],p[i]])
                else:
                    luminal_dataPoints.append([p[0], p[i]])
 
        #shuffle data 
        basal_dataPoints, luminal_dataPoints = shuffle(basal_dataPoints, luminal_dataPoints)
        self.trainBasal = basal_dataPoints[0:len(basal_dataPoints)*.7]
        self.testBasal = basal_dataPoints[len(basal_dataPoints)*.7:]
        self.trainLuminal = luminal_dataPoints[0:len(luminal_dataPoints)*.7]
        self.testLuminal = luminal_dataPoints[len(luminal_dataPoints)*.7:]

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


    ########### Luminal Functions #############
    #leave out test set
    def get_luminal_training_data(self):
        return self.trainLuminal


    def get_luminal_testing_data(self):
        pass

    def get_y_luminal_data(self):
        pass
    
