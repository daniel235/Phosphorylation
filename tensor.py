import tensorflow as tf
from sklearn.cluster import KMeans
from sklearn import svm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D


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
        LuminalX = []
        BasalX = []
        ly = []
        by = []

        #start y data
        for i in range(len(pos_of_found_proteins)):
            index = pos_of_found_proteins[i]
            obs = self.objects[index]
            for k in range(2):
                for j in range(len(obs.sites)):
                    #format (protein expression / site expression / luminal or basal)
                    if k == 0:
                        x.append([obs.get_lExpression(), obs.sites[j].get_lExpression()])
                        y.append(obs.sites[j].name)
                        x_counter += 1
                    else:
                        x.append([obs.get_bExpression(), obs.sites[j].get_bExpression()])
                        y.append(obs.sites[j].name)
                        x_counter += 1


        xr, ys = self.prepare_y_data(y)


        #updating x data set with protein sites only found in the kinase substrate data
        #need to grab even mods to luminal
        for i in range(len(xr)):
            if i % 2 == 0:
                LuminalX.append(x[xr[i]])
                ly.append(ys[i])
            else:
                BasalX.append(x[xr[i]])
                by.append(ys[i])


        type_ys = []
        newX.append(LuminalX)
        newX.append(BasalX)
        type_ys.append(ly)
        type_ys.append(by)


        return newX, type_ys, ys



    def split_data(self):
        print(self.data[0][0])
        x = np.array([self.data[0], self.data[2]])
        label = np.array([self.data[1], self.data[3]])
        #todo prepare data set returns site data and expression
        x, type_y, y = self.prepare_x_data(x, label)
        #luminalx / basal x
        lx = x[0]
        bx = x[1]
        ly = type_y[0]
        by = type_y[1]

        #cross validation
        kfold = KFold(3, True, 1)
        Luminal_three_fold = []
        Basal_three_fold = []
        xtestData = []
        xtrainData = []
        ytrainData = []
        ytestData = []

        #one hot y data as sparse matrix
        labels, uniques = pd.factorize(y)


        labels = labels.tolist()
        uniques = uniques.tolist()

        #creating train and test data to pass to svm
        #Luminal Data
        for train, test in kfold.split(lx):
            for i in train:
                xtrainData.append(lx[i])
                ytrainData.append(labels[labels.index(uniques.index(ly[i]))])

            for j in test:
                xtestData.append(lx[j])
                ytestData.append(labels[labels.index(uniques.index(ly[j]))])


            Luminal_three_fold.append([xtrainData, ytrainData, xtestData, ytestData])
            xtrainData = []
            xtestData = []
            ytrainData = []
            ytestData = []

        #Basal Data
        for train, test in kfold.split(bx):
            for i in train:
                xtrainData.append(bx[i])
                ytrainData.append(labels[labels.index(uniques.index(by[i]))])

            for j in test:
                xtestData.append(bx[j])
                ytestData.append(labels[labels.index(uniques.index(by[j]))])


            Basal_three_fold.append([xtrainData, ytrainData, xtestData, ytestData])
            xtrainData = []
            xtestData = []
            ytrainData = []
            ytestData = []

        return Luminal_three_fold, Basal_three_fold


    def plot_data(self, x, y):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        xs = []
        ys = []
        zs = y

        for i in range(len(x)):
            xs.append(x[i][0])
            ys.append(x[i][1])

        ax.set_xlabel("Protein Expression")
        ax.set_ylabel("Site Expression")

        ax.scatter(xs, ys, zs)
        plt.show()

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
        luminal_data, basal_data = self.split_data()

        accuracy = 0
        #todo separate luminal and basal svm
        #todo mod numbers to get different data sets
        #luminal svm
        clf = SVC(kernel="poly", gamma='scale', verbose=1)
        for i in range(len(luminal_data)):
            clf.fit(luminal_data[i][0], luminal_data[i][1])

        for i in range(len(luminal_data)):
            for j in range(len(luminal_data[i][2])):
                if clf.predict([luminal_data[i][2][j]]) == luminal_data[i][3][j]:
                    accuracy += 1


            print("luminal accuracy ", i, (accuracy / len(luminal_data[i][2])))
            accuracy = 0


        self.plot_data(luminal_data[0][0], luminal_data[0][1])

        #basal svm
        for i in range(len(basal_data)):
            clf.fit(basal_data[i][0], basal_data[i][1])

        for i in range(len(basal_data)):
            for j in range(len(basal_data[i][2])):
                if clf.predict([basal_data[i][2][j]]) == basal_data[i][3][j]:
                    accuracy += 1


            print("basal accuracy ", i, (accuracy / len(basal_data[i][2])))
            accuracy = 0


        self.plot_data(luminal_data[0][0], luminal_data[0][1])



        #get accuracy
        accuracy = 0
        correct = 0
        wrong = 0


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