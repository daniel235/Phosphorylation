import pandas as pd
import numpy as np
class SubstrateCounts:
    def __init__(self, filename=None):
        self.filename = filename


    def get_substrate_counts_in_data(self):
        #import file
        substrates = np.array(pd.read_csv(self.filename, delimiter="\t"))
        substrates = substrates[:,1]
        print(substrates)


s = SubstrateCounts("./data/BC_1.txt")
s.get_substrate_counts_in_data()