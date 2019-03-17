import numpy as np
import pandas as pd
import tensor
import pipe_line as pipe
import protein_interaction_predictor as protein_interaction

#grab data
#package data together

pipe_object = pipe.Pipe_line()
data = pipe_object.get_data()

proteins = []
d = np.array([data[2], data[0]])


#array of protein objects
protein_objects = pipe_object.find_matching_data(d)

#protein interaction network
pInteract = protein_interaction.protein_interaction_net(protein_objects)
#pInteract.network()




#start network call
model = tensor.Network(data, protein_objects, pipe_object)
#model.cluster_network()

model.split_data()
model.regression_network(epochs=10000)



