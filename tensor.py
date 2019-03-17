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


    def softmax_classification_network(self):
        X = tf.placeholder(dtype="tf.float32", shape=(None, 2))
        b = tf.placeholder(dtype="tf.float32", shape=(None, 1))
        Y = tf.placeholder(dtype="tf.int32", shape=(1,1))

        W1 = tf.placeholder(shape=(2, 1), type="tf.float32")

        logit1 = tf.matmul(X * W1) + b
        layer1 = tf.nn.sigmoid(logit1)

        W2 = tf.placeholder(dtype="tf.float32", shape=(1,1))
        b2 = tf.placeholder(dtype="tf.float32", shape=(1,1))

        logit2 = tf.matmul(layer1 * W2) + b2
        layer2 = tf.nn.relu(logit2)

        #layer 3
        W3 = tf.placeholder(dtyype="tf.float32", shape=(1,1))
        b3 = tf.placeholder(dtype="tf.float32", shape=(1,1))

        logit3 = tf.matmul(layer2 * W3) + b3
        layer3 = tf.nn.elu(logit3)

        #output layer
        W4 = tf.placeholder(dtype="tf.float32", shape=(1, 52))
        b4 = tf.placeholder(dtype="tf.float32", shape=(1, 52))


        logit4 = tf.matmul(layer3 * W4) + b4
        layer4 = tf.nn.softmax(logit4)

        #grab data


        loss = tf.reduce_mean(Y - tf.argmax(layer4))
        print(loss)
        #multi class -> softmax = tf.nn.softmax_cross_entropy_with_logits(layer2)
        train = tf.train.GradientDescentOptimizer(learning_rate=.01)
        trainer = train.minimize(loss)

        with tf.Session() as sess:
            tf.global_variables_initializer()
            #sess.run(trainer, feed_dict={X: xData, Y: yData})

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

        #get accuracy
        accuracy = 0
        correct = 0
        wrong = 0

    def pipe_line_svm(self, data):
        clf = svm.SVC(kernel="Linear")
        pipeline = Pipeline(['svc', clf])
        pipeline.fit(data[0], data[1])


    def train_network(self, layer):
        #parameters
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
        cost = tf.reduce_sum(layer)


    def convert_data_for_softmax_network(self):
        pass
        #shape for x -> (None, 2)
        #shape for y -> (1, 1)


    #function for getting impact of different parameters
    def parameter_effects(self):
        #regression network
        pass
        #least squares approach

    def regression_network(self, learningRate=None, epochs=None):
        xLen = len(self.luminal_trainX)
        print(self.luminal_trainX[0])
        xdata = np.array(self.luminal_trainX)
        #xdata = np.reshape(xdata, (2,103))
        #xdata = xdata.transpose()

        print(xdata.shape)


        x = tf.placeholder(dtype=tf.float32, name="input", shape=(1,2))
        y= tf.placeholder(dtype=tf.float32, name="output")


        #coefficients
        W = tf.Variable(initial_value=np.random.uniform(size=(2,1)), name="slope", dtype=tf.float32)
        b = tf.Variable(initial_value=np.random.uniform(size=(1,1)), name="bias", dtype=tf.float32)

        learning_rate = 0.01
        if learningRate != None:
            learning_rate = learningRate


        training_epochs = 2000
        if epochs != None:
            training_epochs = epochs



        #multiplication should be -> (132, 2) * (2, 2) / reshape x to (2, 132) weight should be (132, 1)  outcome would be (2, 1)
        y_pred = tf.add(tf.matmul(x, W), b)


        self.luminal_trainY = np.array(self.luminal_trainY)
        self.luminal_trainY = self.luminal_trainY.reshape((103, 1))

        #self.luminal_trainY.astype(float)

        #mean squared error
        cost = tf.reduce_sum(tf.pow(y_pred - y, 2)) / (2 * xLen)


        #gradient descent optimizer
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=learning_rate).minimize(cost)

        init = tf.global_variables_initializer()

        accuracy = 0

        print(len(xdata))

        with tf.Session() as sess:
            sess.run(init)

            for epoch in range(training_epochs):
                for _x, _y in zip(xdata, self.luminal_trainY):
                    _x = np.reshape(_x, (1, 2))
                    sess.run(optimizer, feed_dict={x: _x, y: _y})

                #sess.run(optimizer, feed_dict={x: xdata, y: self.luminal_trainY})

                if (epoch + 1) % 50 == 0:
                    for _x, _y in zip(xdata, self.luminal_trainY):
                        _x = np.reshape(_x, (1, 2))
                        c = sess.run(cost, feed_dict={x: _x, y: _y})
                        accuracy += c

                        #c = sess.run(cost, feed_dict={x: xdata, y: self.luminal_trainY})
                        print("epoch", epoch, " cost ", c, " W ", sess.run(W), " b ", sess.run(b))

                    accuracy = accuracy / len(xdata)
                    print("accuracy ", (1 - accuracy))


            #training_cost = sess.run(cost, feed_dict={x: xdata, y: self.luminal_trainY})
            Weight = sess.run(W)
            bias = sess.run(b)

            #print(training_cost, Weight, bias)


            #############printing all tests ###############
            print(bias)


            #plot data
            '''print(type(Weight))
            Weight = int(Weight)
            bias = float(bias)
            guess = []
            for i in range(len(self.luminal_trainX)):
                guess.append((Weight * self.luminal_trainX[i]) + bias)

            plt.plot(self.luminal_trainX, self.luminal_trainY, 'ro')
            plt.plot(self.luminal_trainX, guess)
            plt.show()'''
