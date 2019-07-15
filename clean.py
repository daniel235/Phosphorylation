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
        self.data = np.delete(self.data, nth_column, 1)
       