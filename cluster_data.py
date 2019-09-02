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

import clean
import stats


class ClusterData:
    '''This class creates the prepared data for the knn and also hierarchy clustering
    functions.  

    Reads in Excel File

    All functions are separated by underscores

    
    Returns:
        prepared data for clustering algorithms as class properties -> (self.trainBasal)
    '''

    def __init__(self, kinaseSubstrateFile):
        self.pfile = None
        self.fileName = None
        self.kfile = os.path.join(kinaseSubstrateFile)
        self.strongKinase = []
        self.weakKinase = []
        self.kinaseCounts = []
        self.kinaseData = np.array(pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t"))
        self.phosphorylationData = np.array(pd.read_csv("./data/phosphorylation_data.txt", delimiter="\t"))
        self.CancerData = None
        #self.breastCancerData = xlrd.open_workbook("./data/BreastCancerData.xlsx", encoding_override="cp1252").sheet_by_name("data")
        self.phosphositePlusKinaseData = np.array(pd.read_csv("./data/KSA_human.txt", delim_whitespace=True))
        self.unique_kinases = None
        self.colNames = None
        self.phosDataOrdered = True
        self.stats = stats.Statistics()
        self.clean_data()
        

    def replace_with_average(self):
        delete_rows = []
        #for every element in array with na replace with average
        #first get average of row
        for i in range(len(self.CancerData[:,1])):
            indexes = []
            average = 0
            non_empty = 0

            for j in range(2, len(self.CancerData[i])):
                if self.fileName == 'colorectal_cancer.xlsx':
                    if type(self.CancerData[i,j]) != float:
                        indexes.append(j)

                elif np.isnan(self.CancerData[i,j]):
                    indexes.append(j)

                else:
                    average += self.CancerData[i,j]
                    non_empty += 1

            average = average / max(1, non_empty)

            #delete row
            if float(len(indexes) / len(self.CancerData[i])) > .5:
                #replace k with k+1
                delete_rows.append(i)
            else:
                for k in indexes:
                    self.CancerData[i,k] = average

        self.CancerData = np.delete(self.CancerData, delete_rows, 0)


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
        self.fileName = input("What file do you want to use?")
        self.pfile = os.path.join("./data/", self.fileName)
        sheet_name = input("What is your sheet name for phosphorylation data?")
        inputs = input("Is your phosphorylation data ordered(yes/no)?")
        if inputs == "yes":
            print("yes ordered")
            self.phosDataOrdered = True
        else:
            self.phosDataOrdered = False

        inputs = input("Trailing letters on genes?(yes/no)")
        if inputs == "yes":
            trailing = True
        else:
            trailing = False

        df = clean.cleanMatrix(self.pfile, sheet_name)

        #set stats (len of psites, len of unique kinases, len of tumor samples)
        self.stats.set_table(len(df.data[:,0]), len(self.unique_kinases), len(df.data[1]))
        self.stats.plotTable()

        input_column = input("Which column(s) is your psite in?(separate by space if more than 1)")
        index = []
        for i in range(len(input_column)):
            if input_column[i] != ' ':
                index.append(int(input_column[i]))
        

        df.set_gene_site_column(index, trailing)
        index = []
        current_int = []
        counter = 0
        omit_column = input("what column(s) do i omit?(separate by space)")
        for i in range(len(omit_column)):
            if omit_column[i] != ' ':
                current_int.append(omit_column[i])
            
            else:
                if len(current_int) > 1:
                    for j in range(len(current_int)):
                        if j == 0:
                            current_int[j] = 10 * int(current_int[j])

                        index.append(int(current_int[j]) + current_int[0])
                else:
                    index.append(current_int[0])

                current_int = []

        df.omit_columns(index)
        df.column_check_strings()
        df.clean_rows()
    

        ''' self.CancerData = np.array(pd.read_excel(self.pfile, sheet_name="data", dtype=object))
        #join first two columns
        for i in range(len(self.CancerData[:,0])):
            self.CancerData[i,0] = str(self.CancerData[i,0]) + '-' +  str(self.CancerData[i, 1])
            #strip last letter
            self.CancerData[i,0] = (self.CancerData[i,0])[0:-2]
        '''
        self.CancerData = df.data
        #fix Na's here
        print(self.CancerData[4])
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

        for i in range(len(self.phosphositePlusKinaseData)): 
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
                mid = int(len(self.CancerData) / 2)
                current = mid
                start = True
                prevMid = 0
                substrate = self.phosphositePlusKinaseData[i][1].upper()
                
                #find substrate in phosphorylation data
                #binary search the data
                if PhosDataOrdered:
                    mins = 0
                    while(mid > 0 and prevMid != mid):
                        #if not found
                        if prevMid == mid:
                            break

                        #substrate at current index
                        index = self.CancerData[int(current), 0]
                        if type(index) == float:
                            break

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
                            count += 1
                            #add substrate to names and add it's data row to matrix
                            substrate_names.append(substrate)
                            data = []
                            for i in range(2, len(self.CancerData[current])):
                                data.append(self.CancerData[current][i])

                            substrate_matrix.append(data)
                            break

                #brute force data
                else:
                    for i in range(len(self.CancerData[:,0])):
                        index = self.CancerData[i, 0]
                        if type(index) == float:
                            break

                        #didn't find substrate
                        if substrate == index:
                            substrate_names.append(substrate)
                            data = []
                            for j in range(2, len(self.CancerData[i])):
                                data.append(self.CancerData[i][j])

                            
                            substrate_matrix.append(data)
                            break

                        
            #this should stop loop once kinase name is passed (since its alphabetical)
            elif(start and kinaseFileOrdered and self.phosphositePlusKinaseData[i][0] != kinase):
                break

        return substrate_names, substrate_matrix

                
    #returns kinase matrix dictionary // Kinase is for testing purposes
    def get_kinase_substrate_matrixes(self, threshold, kinase=None):
        kinase_matrixes = {}
        substrates = {}
        names = []
        data = []
        

        #kinases = list(set(self.kinaseData[:,0]))
        #new kinase data
        kinases = np.array(self.unique_kinases)
        #todo hierarcharl clustering
        
        with open("ksa2.txt", 'w+') as f:
            f.write(str(self.pfile))
            for i in range(len(kinases)):
                count = self.count_substrates(kinases[i], ordered=False)
                substrates = {}
                names = []
                if count >= threshold: 
                    names, data = self.grab_substrates(kinases[i], False, PhosDataOrdered=self.phosDataOrdered)
                    if len(names) > threshold:
                        for j in range(len(names)):
                            substrates[names[j]] = data[j]
                            kinase_matrixes[kinases[i]] = substrates
                            

                        f.write(F'{kinases[i]}  {list(substrates.keys())}' + "\n")

        return kinase_matrixes



