import numpy as np
import pandas as pd
import tensor
import pipe_line as pipe

#grab data
#package data together

data = pipe.Pipe_line()
data = data.get_data()

#start network call
model = tensor.Network(data)
model.cluster_network()



