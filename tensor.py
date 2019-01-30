import tensorflow as tf
from sklearn.cluster import KMeans
from sklearn import svm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split


class Network:
    def __init__(self, data, objects, pipe):
        self.data = np.array(data)
        self.objects = np.array(objects)
        self.pipe = pipe


    def prepare_y_data(self, y):
        #get names of kinase data
        y_kinase = []
        x = []
        ys = []
        kinase = self.data[6]
        pos = 0
        substrate = kinase["Substrate"].tolist()
        for i in range(len(y)):
            try:
                if y[i] in substrate:
                    pos = substrate.index(y[i].upper())
                    x.append(i)
                    ys.append(kinase["Kinase"][pos])

            except:
                pass


        return x, ys


    #todo preparing for svm
    def prepare_x_data(self, data, label):
        #get object names in a list  #todo should be done in [find_matching_data]
        names = []
        pos = 0
        x_counter = 0
        x = []
        pos_counter = 0
        pos_of_found_proteins = []

        obs = self.objects.tolist()
        for k in range(len(obs)):
            names.append(obs[k].name)

        #get dimension of data
        lengths = []
        for i in range(len(data)):
            lengths.append(len(data[i]))

        for i in range(len(lengths)):
            for j in range(lengths[i]):
                if i == 0:
                    #todo push phosphosites to protein object
                    #strip last hyphen and numbers
                    protein = data[i][j][0]
                    protein = self.pipe.strip_sites(protein)

                    try:
                        pos = names.index("'" + protein + "'")
                        if pos not in pos_of_found_proteins:
                            pos_of_found_proteins.append(pos)
                        if label[i][j] == 'Luminal':
                            self.objects[pos].add_sites(data[i][j][0], lexp=data[i][j][1])

                        else:
                            self.objects[pos].add_sites(data[i][j][0], bexpression=data[i][j][1])

                    except:
                        pass
                #set x's here
                else:
                    #todo push expression to protein object
                    try:
                        protein = data[i][j][0]
                        pos = names.index(protein)
                        if label[i][j] == 'Luminal':
                            self.objects[pos].lExpressionSum += data[i][j][1]
                            self.objects[pos].l_expression_count += 1

                        else:
                            self.objects[pos].bExpressionSum += data[i][j][1]
                            self.objects[pos].b_expression_count += 1

                    except:
                        pass

        #corresponding protein site for x data
        kinase = self.data[6]

        y = []
        newX = []
        #start y data
        for i in range(len(pos_of_found_proteins)):
            index = pos_of_found_proteins[i]
            obs = self.objects[index]
            for k in range(2):
                for j in range(len(obs.sites)):
                    #format (protein expression / site expression / luminal or basal)
                    if k == 0:
                        x.append([obs.get_lExpression(), obs.sites[j].get_lExpression(), 0])
                        y.append(obs.sites[j].name)
                        x_counter += 1
                    else:
                        x.append([obs.get_bExpression(), obs.sites[j].get_bExpression(), 1])
                        y.append(obs.sites[j].name)
                        x_counter += 1


        xr, ys = self.prepare_y_data(y)

        for i in range(len(xr)):
            newX.append(x[xr[i]])

        return newX, ys



    def split_data(self, randomSeed):
        print(self.data[0][0])
        x = np.array([self.data[0], self.data[2]])
        label = np.array([self.data[1], self.data[3]])
        #todo prepare data set returns site data and expression
        x, y = self.prepare_x_data(x, label)
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=randomSeed)
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
        x_train, x_test, y_train, y_test = self.split_data(0)
        x = [[1], [2], [3]]
        y = [1, 2, 3]

        clf = svm.SVC(kernel="poly", gamma='scale', verbose=1)
        clf.fit(x_train, y_train)
        print(x_test[0], y_test[0])
        print("result ", clf.predict([x_test[0]]))

        #get accuracy
        accuracy = 0
        correct = 0
        wrong = 0


        for j in range(2):
            if j == 1:
                x_train, x_test, y_train, y_test = self.split_data(1)
                clf.fit(x_train, y_train)

            for i in range(len(x_test)):
                if clf.predict([x_test[i]]) == y_test[i]:
                    accuracy += 1
                    correct += 1
                    print("correct")

                else:
                    wrong += 1

            print("correct ", correct, " wrong ", wrong)


        accuracy = accuracy / len(x_test)
        print("accuracy ", accuracy)

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