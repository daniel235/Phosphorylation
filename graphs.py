#this file is to visualize all my matrixes graphs and list of data etc.

def getKinaseFeatures(ks):
    with open("kfeats.txt", 'w+') as f:
        for k in ks:
            f.write(k)