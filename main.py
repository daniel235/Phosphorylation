import tensorflow as tf
import pandas as pd


kinaseData = pd.read_csv('./data/Kinase_Substrates.txt', delimiter="\t")

phosphorylation = pd.read_csv('./data/phosphorylation_data.txt', delimiter="\t")
print(phosphorylation)
