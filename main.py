import pca
import numpy as numpy
import pandas as pd
import clean

'''Get PCA of entire breast cancer dataset'''
#import file

data = clean.cleanMatrix('./data/BreastCancerData.xlsx', 'data')
data.set_gene_site_column([0, 1], True)
data.omit_columns([16, 17, 18])
data.clean_rows()
data = data.data

#?get pca
#remove column name

print(type(data[:,1].tolist()))
vec, u, s, vt = pca.getFeatureVector(data[:,1:], 1)
#print(vec)
