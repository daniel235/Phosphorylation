import tensorflow as tf
from sklearn.cluster import KMeans
from sklearn import svm
import numpy as np


class Network:
    def __init__(self, data):
        self.data = data

    def regression_network(self):
        X = tf.placeholder(dtype="tf.float32", shape=(1, 2))
        b = tf.placeholder(dtype="tf.float32", shape=(1,1))
        Y = tf.placeholder(dtype="tf.int32", shape=(1,1))

        W1 = tf.placeholder(shape=(2, 1), type="tf.float32")

        logit1 = tf.matmul(X * W1) + b
        layer1 = tf.nn.relu(logit1)

        W2 = tf.placeholder(dtype="tf.float32", shape=(1,1))
        b2 = tf.placeholder(dtype="tf.float32", shape=(1,1))


        logit2 = tf.matmul(X * W2) + b2
        layer2 = tf.nn.relu(logit2)

        softmax = tf.nn.softmax_cross_entropy_with_logits(layer2)


    def cluster_network(self):
        #classes
        k = 10
        #how to get all data points on graph
        


        kmeans = KMeans(n_clusters=k)
        kmeans.fit()



    def train_network(self, layer):
        #parameters
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
        cost = tf.reduce_sum(layer)