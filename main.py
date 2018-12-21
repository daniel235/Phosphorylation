import numpy as np
import pandas as pd
import tensor
import pipe_line

#grab data
#package data together
data = pipe_line.Pipe_line.get_data()

#start network call
model = tensor.Network(data)
#model.cluster_network()



