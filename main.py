import numpy as np
import pandas as pd
import network
import pipe_line as pipe
import protein_interaction_predictor as protein_interaction
import svm
import cluster_data as cd
import kmeans as k

#grab data
#package data together

pipe_object = pipe.Pipe_line()
data = pipe_object.get_data()

proteins = []
d = np.array([data[2], data[0]])

#array of protein objects
protein_objects = pipe_object.find_matching_data(d)

pipe_object.grab_substrates('EIF2AK1')


startCluster = cd.ClusterData()
startCluster.set_arbitrary_kinase_class(3)
#startCluster.get_basal_bicor_correlation_matrix()
#k.send_file()
#kmeans = k.Kmeans_cluster()
#kmeans.run_kmeans()

'''
#protein interaction network
pInteract = protein_interaction.protein_interaction_net(protein_objects)
#pInteract.network()

#start network call
model = tensor.Network(data, protein_objects, pipe_object)
#model.cluster_network()

#call regression network to estimate paramters
model.regression_network()'''


