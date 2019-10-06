import pandas as pd
import numpy as np
import pickle


#read unique kinases
filename = "./data/pickles/uniqueKinases"
with open(filename, 'rb+') as f:
    uniqueKinases = pickle.load(f)

#create matrix
uniqueKinases = np.transpose(uniqueKinases)
interaction_matrix = pd.DataFrame(columns=uniqueKinases)
print(interaction_matrix)

#get correlation values
