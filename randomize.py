import pandas as pd
import numpy as np
import clean
#import data

cancer_data = np.array(pd.read_excel("./data/BreastCancerData.xlsx", sheet_name="data", dtype=object))
print(cancer_data)

#clean data
cleanData = clean.cleanMatrix()
cleanData.data = cancer_data
cleanData.clean_rows()
cleanData.set_gene_site_column([0,1], True)


