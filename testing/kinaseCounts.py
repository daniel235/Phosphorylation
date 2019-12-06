import numpy as np
import pandas as pd


class Kinase_Count_Test:
    '''
    Used to make sure we have accurate count of kinases in the file
    '''
    def __init__(self, file):
        self.count = None
        self.fileName = file


    def test_unique_kinases(self, uniqueKinases):
        set1 = uniqueKinases
        sub_data = pd.read_csv("./data/KSA_human.txt", delim_whitespace=True)
        sub_data = np.array(sub_data)
        for i in range(len(sub_data)):
            sub_data[i][1] = sub_data[i][1] + "-" + str(sub_data[i][2])

        sub_data = np.delete(sub_data, 2, 1)

        set2 = set(list(sub_data[:,1]))

        finalCount = set1.intersection(set2)
        return finalCount

    
    def 