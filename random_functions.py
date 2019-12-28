import numpy as np
import pickle

def get_sig_scores(obj, cg):
    obj = np.array(obj)
    

    ##remove zeros from low list
    for k in range(len(obj)):
        if obj[k] == 0:
            #add 5 to it
            obj[k] = 5
    

    o = obj.argsort()
    indexes = []
    lowScores = []
    for p in range(5):
        print(obj[o[p]])
        lowScores.append(obj[o[p]])
        indexes.append(o[p])
        

    print(indexes)
    sigNodes = []
    for j in range(len(indexes)):
        index1 = int(math.floor(indexes[j] / 13))
        index2 = int(math.floor(indexes[j] % 13)) 
        #!todo not going to find it
        sigNodes.append([cg[1][index1].name, cg[0][index2].name])

    #write sig nodes to pickle
    with open("./data/pickles/randomSignificantNodes", 'wb+') as f:
        pickle.dump(sigNodes, f)

    print(sigNodes)

def distro_scores(scores):
    sns.distplot(scores)