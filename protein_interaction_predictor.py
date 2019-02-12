import numpy as np
import pandas as pd

class protein_interaction_net():
    def __init__(self, objects):
        self.protein_objects = objects
        self.protein_interaction_data = pd.read_csv('./data/Protein_Protein_Interaction.txt', delimiter="\t")



