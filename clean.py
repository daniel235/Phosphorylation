import pandas as pd
import numpy as np
import math


class cleanMatrix:
    '''
        this class prepares the data to remove unwanted columns
        prepares the gene site column as one and without tumor type

        Formula 
        1.) omit_columns()
        2.) set_gene_site_column()
        3.) clean_data()

    '''
    def __init__(self, dataFile=None, sheetName=None):
        self.matrix = []
        self.data = None
        if dataFile != None:
            try:
                self.data = pd.read_excel(dataFile, sheet_name=sheetName, dtype=object)
            except:
                self.data = pd.read_csv(dataFile, delim_whitespace=True)
                self.data = np.array(self.data)
                for i in range(len(self.data)):
                    #remove quote marks
                    self.data[i][0] = self.data[i][0].replace("'", '')

                return

        self.data = np.array(self.data)
        

    #remove columns from dataset
    def omit_columns(self, nth_columns):
        '''Delete columns from data
            Example: omit_columns(nth_columns=[1,2,3])
        '''
        self.data = np.delete(self.data, nth_columns, 1)


    #check if columns contain strings and convert to floats
    def column_check_strings(self):
        '''Changes excel string values to floats'''
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
        '''Set the main psite column with complete protein and site name
            nth_columns needs to be in array format nth_columns=[1,2]
            optional: trailing letter to remove T/F
            Also deletes first row to replace
        '''
        #combine columns and if trailing letter(s,t) remove it.
        #first column should be site name
        if len(nth_columns) == 1:
            for i in range(len(self.data[:,0])):    
                self.data[i,0] = self.data[i, nth_columns[0]]

                if trailing_letter:
                    self.data[i,0] = str(self.data[i,0])[:-1]
                    #also check for semicolons
                    #self.check_for_multiples(self.data[i,0])
                    

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

            elif string[i] == ';' or string[i] == ' ' or string[i] == '\t':
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
        '''Deletes rows with more than 50 percent of data missing '''
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
