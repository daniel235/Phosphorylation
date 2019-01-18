import tensorflow as tf
from sklearn.cluster import KMeans
from sklearn import svm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

#todo onehot data


class Network:
    def __init__(self, data, objects, pipe):
        self.data = np.array(data)
        self.objects = np.array(objects)
        self.pipe = pipe


    def prepare_data(self, data):
        #get object names in a list  #todo should be done in [find_matching_data]
        names = []
        obs = self.objects.tolist()
        for k in range(len(obs)):
            names.append(obs[k].name)

        #get dimension of data
        lengths = []
        for i in range(len(data)):
            lengths.append(len(data[i]))
            print(len(data[i]))

        for i in range(len(lengths)):
            for j in range(lengths[i]):
                if i == 0:
                    #todo push phosphosites to protein object
                    #strip last hyphen and numbers
                    protein = data[i][j][0]
                    protein = self.pipe.strip_sites(protein)

                    try:
                        pos = names.index("'" + protein + "'")
                        self.objects[pos].add_sites(data[i][j][0])
                        print("sites", self.objects[j].sites)
                    except:
                        pass

                else:
                    #todo push expression to protein object
                    pass


    def split_data(self):
        print(self.data[0][0])
        x = np.array([self.data[0], self.data[2]])
        y = np.array([self.data[1], self.data[3]])
        self.prepare_data(x)
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=0)
        print(X_train)
        return X_train, X_test, y_train, y_test

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

        loss = tf.reduce_mean(Y - logit2)
        #multi class -> softmax = tf.nn.softmax_cross_entropy_with_logits(layer2)
        train = tf.train.GradientDescentOptimizer(learning_rate=.01)
        train.minimize(loss)


    def cluster_network(self):
        #classes
        #how to get all data points on graph
        x, x2, y, y2 = self.split_data()
        print(x)
        #plt.plot(x, y, 'o')
        #plt.show()

        clf = svm.SVC(kernel="poly", gamma='scale', verbose=1)
        clf.fit(x, y)

        #print("got to fit ")
        #print(clf.fit(x, y))


    def train_network(self, layer):
        #parameters
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
        cost = tf.reduce_sum(layer)


    def one_hot(self, x, y):
        taken_names = []
        y_names = []
        newX = []
        newY = []
        counter = 0
        pos = 0
        for i in range(len(x)):
            if x[i] not in taken_names:
                taken_names.append(x[i])
                newX.append(counter)
            else:
                pos = taken_names.index(x[i])
                newX.append(pos)


        counter = 0
        for j in range(len(y)):
            if y[j] not in y_names:
                y_names.append(y[j])
                newY.append(counter)