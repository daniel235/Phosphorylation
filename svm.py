import tensorflow as tf 
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.svm import LinearSVC, NuSVC
import numpy as np
import pipe_line as pipe 


class SVM:
    #?class to classify phosphosites to kinases
    #?using svm classifier
    def __init__(self):
        pass

    def run_svm(self, test=False):
        svm = NuSVC(decision_function_shape="ovr", verbose=True)
        #?test svm
        if(test):
            #create sample data (3 classes)
            x = [[1,2],[1,3],[1,4],[1,1],[5,1],[5,2],[5,3],[5,4],[7,5],[7,9],[6,8],[7,4]]
            y = [1,1,1,1,2,2,2,2,3,3,3,3]
            svm.fit(x,y)
            pred = svm.predict([[8,2]])
            print(pred)
        else:
            #load data
            pass

    
    #graphing phosphorylation data by different colors
    def graph_data_by_kinases(self):
        pass