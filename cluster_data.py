import numpy as np
import pipe_line as pipe


class clusterData:
    def __init__(self):
        self.basal_data = None
        self.pipe_object = pipe.Pipe_line()

        
    #grab training data from all databases
    def get_basal_training_data(self):
        phosphosite_data = np.array(self.pipe_object.phosphorylation)
        dataPoints = []
        
        for count in range(len(phosphosite_data)):
            p = phosphosite_data[count]
            for i in range(1, 18):
                if i <= 9:
                    dataPoints.append([p[0],p[i]])

        return dataPoints

    #leave out test set
    def get_basal_testing_data(self):
        pass


    def get_luminal_training_data(self): 
        phosphosite_data = np.array(self.pipe_object.phosphorylation)
        dataPoints = []

        for count in range(len(phosphosite_data)):
            p = phosphosite_data[count]
            for i in range(1, 18):
                if i > 9:
                    dataPoints.append([p[0],p[i]])

        return dataPoints


    def get_luminal_testing_data(self):
        pass
