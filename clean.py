import pandas as pd
import numpy as np


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
        

    def omit_columns(self, nth_columns):
        self.data = np.delete(self.data, nth_columns, 1)


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
                
            print("first one ", self.data[0])
            return


    def check_for_multiples(self, string):
        #check for semicolon
        #add to array if multiple
        #ARGTF-567;855;320;
        for i in range(len(string)):
            if string[i] == ';':
                #iterate backwords to save word
                pass
        pass

       