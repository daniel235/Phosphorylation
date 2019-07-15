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
        print(self.data)


    def set_gene_site_column(self, nth_columns, trailing_letter=False):
        #combine columns and if trailing letter(s,t) remove it.
        #first column should be site name
        print(nth_columns)
        if len(nth_columns) == 1:
            self.data[:,0] = self.data[:,nth_columns[0]]
            if trailing_letter:
                for i in range(len(self.data[:,0])):
                    self.data[i,0] = str(self.data[i,0])[:-1]

        else:
            self.data[:,0] = str(self.data[:,nth_columns[0]]) + "-" +  str(self.data[:,nth_columns[1]])
            self.data = np.delete(self.data, nth_columns[1:])
            return

        #drop nth column
        self.data = np.delete(self.data, nth_columns[1:], 1)