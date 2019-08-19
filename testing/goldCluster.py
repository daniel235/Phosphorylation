import pandas as pd
import numpy as np

#todo fix file import
class GoldCluster:
    def __init__(self):
        self.filename = None
        self.clusters = None

    def createCluster(self):
        df = pd.read_excel(r"../data/clusterTable.xls", sheet_name="Sheet1")
        print(df)


gc = GoldCluster()
gc.createCluster()