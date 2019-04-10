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
        self.luminal_trainX = None
        self.luminal_trainY = None
        self.basal_trainX = None
        self.basal_trainY = None
        self.split_data()


    def prepare_y_data(self, y):
        #get names of kinase data
        x, ys = [], []

        #create substrate list
        substrate = self.data[6]["Substrate"].tolist()

        #for each substrate site name determine which has a defined kinase
        for i in range(len(y)):
            #if substrate is in kinase file
            try:
                if y[i] in substrate:
                    pos = substrate.index(y[i].upper())
                    x.append(i)
                    ys.append(self.data[6]["Kinase"][pos])

            except:
                pass

        return x, ys


    #todo preparing for svm
    def prepare_x_data(self, data, label):
        #get object names in a list
        names = []
        x_counter = 0
        x = []
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

        one_check = True

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

            #set data here (y is converted to integer)
            if one_check:
                self.luminal_trainX = xtrainData
                self.luminal_trainY = ytrainData
                one_check = False

            xtrainData = []
            xtestData = []
            ytrainData = []
            ytestData = []

        one_check = True

        #Basal Data
        for train, test in kfold.split(bx):
            for i in train:
                xtrainData.append(bx[i])
                ytrainData.append(labels[labels.index(uniques.index(by[i]))])

            for j in test:
                xtestData.append(bx[j])
                ytestData.append(labels[labels.index(uniques.index(by[j]))])


            Basal_three_fold.append([xtrainData, ytrainData, xtestData, ytestData])

            if one_check:
                self.basal_trainX = xtrainData
                self.basal_trainY = ytrainData
                one_check = False

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


    def cluster_network(self, results=False):
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
            #save accuracy
            if results:
                with open("./results/results.txt", "a") as file:
                    file.write("Luminal Accuracy " + str(accuracy / len(luminal_data[i][2])) + "\n")

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

            if results:
                with open("./results/results.txt", "a") as file:
                    file.write("Basal Accuracy " + str(accuracy / len(basal_data[i][2])) + "\n")
                    if i == 2:
                        file.write("\n")
            #save accuracy
            accuracy = 0


        self.plot_data(luminal_data[0][0], luminal_data[0][1])


    def pipe_line_svm(self, data):
        clf = svm.SVC(kernel="Linear")
        pipeline = Pipeline(['svc', clf])
        pipeline.fit(data[0], data[1])


    def regression_network(self, learningRate=None, epochs=None):
        xLen = len(self.luminal_trainX)
        print(self.luminal_trainX[0])
        xdata = np.array(self.luminal_trainX)
        ydata = self.luminal_trainY


        x = tf.placeholder(dtype=tf.float32, name="input", shape=(1,2))
        y= tf.placeholder(dtype=tf.float32, name="output")


        #coefficients
        W = tf.Variable(initial_value=np.random.uniform(size=(2,1)), name="slope", dtype=tf.float32)
        b = tf.Variable(initial_value=np.random.uniform(size=(1,1)), name="bias", dtype=tf.float32)

        learning_rate = 0.01
        if learningRate != None:
            learning_rate = learningRate

        training_epochs = 200
        if epochs != None:
            training_epochs = epochs


        #multiplication should be -> (132, 2) * (2, 2) / reshape x to (2, 132) weight should be (132, 1)  outcome would be (2, 1)
        y_pred = tf.add(tf.matmul(x, W), b)


        self.luminal_trainY = np.array(self.luminal_trainY)
        self.luminal_trainY = self.luminal_trainY.reshape((103, 1))

        #mean squared error
        cost = tf.reduce_sum(tf.pow(y_pred - y, 2)) / (2 * xLen)

        #gradient descent optimizer
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=learning_rate).minimize(cost)

        init = tf.global_variables_initializer()

        accuracy = 0

        ###########  Cross-Fold Data (5x) ##############

        X_train, X_test, y_train, y_test = train_test_split(xdata, ydata, test_size=0.3)
        kf = KFold(n_splits=5)

        ###########  Train Network  ##############


        with tf.Session() as sess:
            sess.run(init)
            newXdata = []
            newYdata = []
            for train, test in kf.split(xdata):
                newXdata.append(xdata[train])
                newYdata.append(ydata[train])
                for epoch in range(training_epochs):
                    for _x, _y in zip(newXdata, newYdata):
                        _x = np.reshape(_x, (1, 2))
                        sess.run(optimizer, feed_dict={x: _x, y: _y})

                    if (epoch + 1) % 50 == 0:
                        for _x, _y in zip(xdata, ydata):
                            _x = np.reshape(_x, (1, 2))
                            c = sess.run(cost, feed_dict={x: _x, y: _y})
                            accuracy += c

                            #c = sess.run(cost, feed_dict={x: xdata, y: self.luminal_trainY})
                            print("epoch", epoch, " cost ", c, " W ", sess.run(W), " b ", sess.run(b))

                        accuracy = accuracy / len(xdata)
                        print("accuracy ", (1 - accuracy))


            Weight = sess.run(W)
            bias = sess.run(b)


            #############printing all tests ###############

            #plot data
            '''guess = []
            for i in range(len(self.luminal_trainX)):
                guess.append((Weight * self.luminal_trainX[i]) + bias)

            plt.plot(self.luminal_trainX, self.luminal_trainY, 'ro')
            plt.plot(self.luminal_trainX, guess)
            plt.show()'''
