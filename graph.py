
        

import pandas as pd

class Graph:
    def __init__(self):
        self.graph = None
        self.kinaseNodes = []
        self.kinaseData = pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t")

    def createNode(self):
        pass

    def add(self):
        pass

    def addEdge(self, node, substrate):
        #add substrate as edge 
        pass
        


    def createGraph(self):
        for i in range(len(self.kinaseData)):
            if self.kinaseData["Kinase"][i] in self.kinaseNodes:
                self.addEdge(self.kinaseNodes[self.kinaseData["Kinase"][i]], self.kinaseData["Substrate"][i])

            else:
                self.createNode()

class Node:
    def __init__(self):
        self.name = ""
        self.edges = []

    def add(self, substrate):
        self.edges.append(substrate)
            
        

