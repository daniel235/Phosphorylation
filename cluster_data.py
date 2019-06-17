import numpy as np
from sklearn.utils import shuffle
import pipe_line as pipe
import math
import pandas as pd
import xlrd 
from scipy.stats import stats
from astropy import stats
from sklearn.decomposition import PCA
import pickle
import os


class ClusterData:
    '''This class creates the prepared data for the knn and also hierarchy clustering
    functions.  

    Reads in Excel File

    All functions are separated by underscores

    returns datasets for basal and luminal types.

    Also creates the y data to test against data.

    Automatically separates data into train and test set

    
    
    Returns:
        prepared data for clustering algorithms as class properties -> (self.trainBasal)
    '''

    def __init__(self, kinaseSubstrateFile, PhosphoDataFile):
        self.pfile = os.path.join(PhosphoDataFile)
        self.kfile = os.path.join(kinaseSubstrateFile)
        self.strongKinase = []
        self.weakKinase = []
        self.kinaseCounts = []
        self.kinaseData = np.array(pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t"))
        self.phosphorylationData = np.array(pd.read_csv("./data/phosphorylation_data.txt", delimiter="\t"))
        self.breastCancerData = None
        #self.breastCancerData = xlrd.open_workbook("./data/BreastCancerData.xlsx", encoding_override="cp1252").sheet_by_name("data")
        self.phosphositePlusKinaseData = np.array(pd.read_csv("./data/KSA_human.txt", delim_whitespace=True))
        self.unique_kinases = None
        self.clean_data()
        
    def replace_with_average(self):
        #for every element in array with na replace with average
        #first get average of row
        for i in range(len(self.breastCancerData)):
            indexes = []
            average = 0
            
            for j in range(2, len(self.breastCancerData[i])):
                if np.isnan(self.breastCancerData[i,j]):
                    indexes.append(j)
                else:
                    average += self.breastCancerData[i,j]

            for k in indexes:
                self.breastCancerData[i,k] = average


    #clean breast cancer data and create 
    #kinase matrix and phosphosite matrix
    def clean_data(self):
        self.unique_kinases = np.array(pd.read_csv(self.kfile, delim_whitespace=True))[:,0]
        unique_kinase_temp = []
        for i in self.unique_kinases:
            if i not in unique_kinase_temp:
                unique_kinase_temp.append(i)
        #set unique kinases
        #strip na's
        #strip columns
        self.unique_kinases = unique_kinase_temp
        
        #print(self.phosphositePlusKinaseData[:,1])
        self.breastCancerData = np.array(pd.read_excel(self.pfile, sheet_name="data", dtype=object))
        #join first two columns
        for i in range(len(self.breastCancerData[:,0])):
            self.breastCancerData[i,0] = str(self.breastCancerData[i,0]) + '-' +  str(self.breastCancerData[i, 1])
            #strip last letter
            self.breastCancerData[i,0] = (self.breastCancerData[i,0])[0:-2]
           
        #fix Na's here
        self.replace_with_average()

        #fix kinase substrates columns
        for i in range(len(self.phosphositePlusKinaseData[:,1])):
            self.phosphositePlusKinaseData[i,1] = str(self.phosphositePlusKinaseData[i,1]) + "-" + str(self.phosphositePlusKinaseData[i,2])
        
        self.phosphositePlusKinaseData = self.phosphositePlusKinaseData[:,0:-1]


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
    
        #separate kinases
        for i in range(len(kinaseDict)):
            if kinaseDict[names[i]] > n:
                self.strongKinase.append(names[i])
            else:
                self.weakKinase.append(names[i])
                
  
    #todo modularize functions
    def count_substrates(self, kinase, ordered=False):
        start = False
        subCount = 0

        for i in range(len(self.kinaseData)): 
            if self.phosphositePlusKinaseData[i][0] == kinase:
                if not start:
                    start = True
                    
                subCount += 1
            #this should stop loop once kinase name is passed (since its alphabetical)
            elif(start and ordered and self.kinaseData[i][0] != kinase):
                break


        return subCount
                
    #?grab substrates of kinase passed in
    def grab_substrates(self, kinase, kinaseFileOrdered=False, PhosDataOrdered=False):
        substrate_matrix = []
        substrate_names = []
        start = False
        count = 0
        for i in range(len(self.phosphositePlusKinaseData)):
            if self.phosphositePlusKinaseData[i][0] == kinase:
                #?logic controllers
                mid = int(len(self.breastCancerData) / 2)
                current = mid
                start = True
                prevMid = 0
                substrate = self.phosphositePlusKinaseData[i][1]
            
                #find substrate in phosphorylation data
                #binary search the data
                if PhosDataOrdered:
                    mins = 0
                    while(mid > 0 and prevMid != mid):
                        #if not found
                        if prevMid == mid:
                            break

                        #substrate at current index
                        index = self.breastCancerData[int(current), 0]
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
                            print("found substrate")
                            count += 1
                            #add substrate to names and add it's data row to matrix
                            substrate_names.append(substrate)
                            data = []
                            for i in range(2, len(self.breastCancerData[current])):
                                data.append(self.breastCancerData[current][i])

                            substrate_matrix.append(data)
                            break
                        
            
            #this should stop loop once kinase name is passed (since its alphabetical)
            elif(start and kinaseFileOrdered and self.phosphositePlusKinaseData[i][0] != kinase):
                break

        print(kinase, " substrate count ", count)
        return substrate_names, substrate_matrix

                
    #returns kinase matrix dictionary
    def get_kinase_substrate_matrixes(self, threshold):
        kinase_matrixes = {}
        substrates = {}
        names = []
        data = []
        substratesLength = []

        #kinases = list(set(self.kinaseData[:,0]))
        #new kinase data
        kinases = self.unique_kinases
        for i in range(len(kinases)):
            count = self.count_substrates(kinases[i], ordered=False)
            substrates = {}
            if count >= threshold:
                names, data = self.grab_substrates(kinases[i], False, PhosDataOrdered=True)
                for i in range(len(names)):
                    if len(names) > threshold:
                        substratesLength.append(len(names[i]))
                        substrates[names[i]] = data[i]
                        kinase_matrixes[kinases[i]] = substrates

        
        return kinase_matrixes, substratesLength


        


