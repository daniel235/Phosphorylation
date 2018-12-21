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
        #how to get all data points on graph
        x = []
        x2 = []
        index = 0
        index2 = 0
        attained = False
        for i in range(len(self.data[3])):
            #find protein index
            for j in range(len(self.data[7])):
                if self.data[3][i][0] == self.data[7][j]:
                    index = j
                    break

            x.append([index, self.data[3][i][1]])

        for i in range(len(self.data[5])):
            for j in range(len(self.data[8])):
                if self.data[5][i][0] == self.data[8][j]:
                    index2 = j
                    break

            print([index2, self.data[5][i][1]])
            x2.append([index2, self.data[5][i][1]])

        x = zip(x, x2)
        x = set(x)
        x = np.array(x)

        y = self.data[4]
        clf = svm.SVC(kernel="poly", gamma='scale', verbose=1)
        print("got to fit ")
        print(clf.fit(x, y))


    def train_network(self, layer):
        #parameters
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
        cost = tf.reduce_sum(layer)