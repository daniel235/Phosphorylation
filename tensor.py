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


    def prepare_data(self, data, label):
        #get object names in a list  #todo should be done in [find_matching_data]
        names = []
        pos = 0
        x_counter = 0
        x = []

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
                        if data[i][j][0] not in self.objects[pos].sites:
                            self.objects[pos].add_sites(data[i][j][0])
                    except:
                        pass
                #set x's here
                else:
                    #todo push expression to protein object
                    try:
                        protein = data[i][j][0]
                        pos = names.index(protein)
                        self.objects[pos].all_expression += data[i][j][1]
                        x[x_counter] = [self.objects[pos].expression, self.objects[pos]]
                        print("expression", self.objects[pos].expression)

                    except:
                        pass


        #start y data





    def split_data(self):
        print(self.data[0][0])
        x = np.array([self.data[0], self.data[2]])
        label = np.array([self.data[1], self.data[3]])
        #todo prepare data set returns site data and expression
        self.prepare_data(x, label)
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
        #x, x2, y, y2 = self.split_data()
        #print(x)
        #plt.plot(x, y, 'o')
        #plt.show()
        self.split_data()
        x = [[1], [2], [3]]
        y = [1, 2, 3]

        clf = svm.SVC(kernel="poly", gamma='scale', verbose=1)
        clf.fit(x, y)
        print("result ", clf.predict([[1]]))

        #print("got to fit ")
        #print(clf.fit(x, y))

    def train_network(self, layer):
        #parameters
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
        cost = tf.reduce_sum(layer)


    #not really one hot
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
                #only increment when new item is added to taken names
                counter += 1
            else:
                #use the position to mark similar x
                pos = taken_names.index(x[i])
                newX.append(pos)

        counter = 0

        for j in range(len(y)):
            if y[j] not in y_names:
                y_names.append(y[j])
                pos = counter
                counter += 1
            else:
                pos = y_names.index(y[j])
                newY.append(pos)

            newY.append(pos)

        return newX, newY