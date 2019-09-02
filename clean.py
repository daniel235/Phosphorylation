import pandas as pd
import numpy as np
import math


class cleanMatrix:
    '''
        this class prepares the data to remove unwanted columns
        prepares the gene site column as one and without tumor type

    '''
    def __init__(self, dataFile, sheetName):
        self.matrix = []
        self.data = pd.read_excel(dataFile, sheet_name=sheetName, dtype=object)
        self.colnames = self.data.columns
        self.data = np.array(self.data)
        

    #remove columns from dataset
    def omit_columns(self, nth_columns):
        self.data = np.delete(self.data, nth_columns, 1)


    #check if columns contain strings and convert to floats
    def column_check_strings(self):
        row = self.data[0]
        change_cols = []
        for i in range(1, len(row)):
            if type(row[i]) != float:
                #convert to float
                change_cols.append(i)


        for i in change_cols:
            for j in range(len(self.data[:,i])):
                try:
                    self.data[j,i] = float(self.data[j,i])
                except:
                    self.data[j,i] = None


    #create the gene site column on the very first column
    def set_gene_site_column(self, nth_columns, trailing_letter=False):
        #combine columns and if trailing letter(s,t) remove it.
        #first column should be site name
        if len(nth_columns) == 1:
            for i in range(len(self.data[:,0])):    
                self.data[i,0] = self.data[i, nth_columns[0]]

                if trailing_letter:
                    self.data[i,0] = str(self.data[i,0])[:-1]
                    #also check for semicolons
                    self.check_for_multiples(self.data[i,0])
                    

            if nth_columns[0] != 0:
                self.data = np.delete(self.data, nth_columns, 1)

        else:
            for i in range(len(self.data[:,0])):
                self.data[i,0] = str(self.data[i,nth_columns[0]]) + "-" +  str(self.data[i,nth_columns[1]])
                if trailing_letter:
                    #self.data[i,0] = str(self.data[i,0])[:-1]
                    while(True):
                        if str(self.data[i,0])[-1] > '9' or str(self.data[i,0])[-1] < '0':
                            self.data[i,0] = str(self.data[i,0])[:-1]
                        else:
                            break

            if nth_columns[0] != 0:
                self.data = np.delete(self.data, nth_columns, 1)

            else:
                self.data = np.delete(self.data, nth_columns[1:], 1)
                

            return


    def check_for_multiples(self, string, index):
        #check for semicolon
        #add to array if multiple
        #ARGTF-567;855;320;
        base = ""
        site = ""
        sites = {}
        start = False
        returnStatus = False
        for i in range(len(string)):
            if string[i] == '-':
                start = True

            elif string[i] == ';':
                #iterate backwords to save word
                sites[base + "-" + site] = index
                returnStatus = True
                site = ""
            
            elif start == True:
                site += string[i]

            else:
                base += string[i]


        if returnStatus == True:
            return sites
        
        else:
            return None
        


    def clean_rows(self):
        delete_rows = []
        add_rows = []
        #for every element in array with na replace with average
        #first get average of row
        #simultaneously get the extra rows and remove added sites
        for i in range(len(self.data[:,1])):
            indexes = []
            average = 0
            non_empty = 0


            for j in range(1, len(self.data[i])):
                try:
                    if np.isnan(self.data[i,j]):
                        indexes.append(j)

                    else:
                        average += self.data[i,j]
                        non_empty += 1

                except:
                    if type(self.data[i,j]) != float:
                        indexes.append(j)
                    
                    else:
                        average += self.data[i,j]
                        non_empty += 1

            average = average / max(1, non_empty)

            #delete row
            if float(len(indexes) / len(self.data[i])) > .5:
                #replace k with k+1
                delete_rows.append(i)
            else:
                for k in indexes:
                    self.data[i,k] = average

        self.data = np.delete(self.data, delete_rows, 0)
       
        #run through all rows and create new rows
        for i in range(len(self.data[:,1])):
            #start here for column site cleanup
            extra_rows = self.check_for_multiples(self.data[i,0], i)
            first = ""
            firstFlag = True

            if extra_rows != None:
                row = self.data[i]
                for name, index in extra_rows.items():
                    if firstFlag == True:
                        first = name
                        firstFlag = False
                    else:
                        #fetch row from data 
                        row[0] = name
                        #add to add_rows
                        add_rows.append(row)

                #set original row name to first extra row
                print("first ", first)
                self.data[i][0] = first


        if len(add_rows) > 1:
            self.data = np.append(self.data, add_rows, axis=0)

