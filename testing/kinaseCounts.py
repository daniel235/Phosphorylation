import numpy as np
import pandas as pd
import sys

#import dev files
sys.path.append("../alias.py")
import alias

class Kinase_Count_Test:
    '''
    Used to make sure we have accurate count of kinases in the file
    '''
    def __init__(self, file=None, alias_object=None):
        self.count = None
        self.fileName = file
        self.alias = alias_object


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

    
    def test_kinase_alias(self, kinase):
        #import alias kinases
        alias_data = pd.read_csv("./data/info_table.csv", delimiter="\t")
        main_alias = list(alias_data["Gene"])
        #upper case kinases
        kinase = kinase.upper()
        for k in range(len(main_alias)):
            main_alias[k] = main_alias[k].upper()

        #first check if kinase exists in alias data
        if kinase not in self.alias.alias_dict.keys():
            return True

        print(main_alias)
        if kinase not in main_alias:
            return False

        return True

