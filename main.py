import numpy as np
import pandas as pd
import tensor
import pipe_line as pipe

#grab data
#package data together

pipe_object = pipe.Pipe_line()
data = pipe_object.get_data()

proteins = []
d = np.array([data[2], data[0]])
pipe_object.find_matching_data(d)



#start network call
#model = tensor.Network(data)
#model.cluster_network()



