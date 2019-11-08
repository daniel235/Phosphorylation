import os
import numpy as np
import pickle


class VisualData:
    def __init__(self):
        self.exitFlag = False
        self.dataObj = None


    def query(self):
        pickles = os.listdir('./data/pickles')
        while(self.exitFlag==False):
            for i in range(len(pickles)):
                print(i+1, " ", pickles[i])
            
            print((len(pickles)+1), " Exit")

            choice = input("What data object do you want?")
            if int(choice) >= len(pickles):
                self.exitFlag = True
                return

            fileName = './data/pickles/' + str(pickles[int(choice) - 1])
            f = open(fileName, 'rb+')
            self.dataObj = pickle.load(f)
            print(self.dataObj)
            self.getType()


    def getType(self):
        if type(self.dataObj) == np.ndarray:
            print(self.dataObj[1,3])
            print(self.dataObj[3,1])
            print(self.dataObj.shape)
        print(type(self.dataObj))
        



v = VisualData()
v.query()