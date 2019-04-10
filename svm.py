import tensorflow as tf 
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.svm import LinearSVC, NuSVC

import pipe_line as pipe 


class SVM:
    #?class to classify phosphosites to kinases
    #?using svm classifier
    def __init__(self):
        self.pipe_object = pipe.Pipe_line()

    #grab training data from all databases
    def get_basal_training_data(self):
        phosphosite_data = np.array(self.pipe_object.phosphorylation)
        dataPoints = []
        
        for count in range(len(phosphosite_data)):
            for i in range(18):


    #leave out test set
    def get_basal_testing_data(self):
        pass

    def get_luminal_training_data(self):
        pass

    def get_luminal_testing_data(self):
        pass


    def run_svm(self):
        svm = NuSVC(decision_function_shape="ovr", verbose=True)
        svm.fit()
    
    #graphing phosphorylation data by different colors
    def graph_data_by_kinases(self):
        pass