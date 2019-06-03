import numpy as np
import pandas as pd
import network
import pipe_line as pipe
import protein_interaction_predictor as protein_interaction
import svm
import cluster_data as cd
import kmeans as k
import graph 

#grab data
#package data together


#graphs = graph.Graph()
#graphs.createGraph()

startCluster = cd.ClusterData()
startCluster.grab_substrates("'CDK1'", True, True)
#startCluster.set_arbitrary_kinase_class(3)
#startCluster.sampleCorrelationMatrix()
#bicor = startCluster.get_basal_bicor_correlation_matrix()
#startCluster.pca(bicor)
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


