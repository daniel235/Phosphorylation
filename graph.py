import pandas as pd

class Graph:
    def __init__(self):
        self.graph = None
        self.kinaseNodes = []
        self.kinaseData = pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t")
        self.edges = []


    def createNode(self, name):
        node = Node(name)
        self.kinaseNodes.append(node)


    def add(self):
        pass


    def addEdge(self, node, substrate):
        #add substrate as edge 
        #check if edge is in use
        if substrate in self.edges:
            index = self.edges.index(substrate)
            #grab neighbor node 
            myNode = self.kinaseNodes[index]

        #no neighbor for node
        else:
            self.edges.append(substrate)
            node.edges.append(substrate)


    def createGraph(self):
        for i in range(len(self.kinaseData)):
            if self.kinaseData["Kinase"][i] in self.kinaseNodes:
                node = self.kinaseNodes.index(self.kinaseData["Kinase"][i])
                self.addEdge(node, self.kinaseData["Substrate"][i])

            else:
                self.createNode(self.kinaseData["Kinase"][i])


class Node:
    def __init__(self, name):
        self.name = name
        self.edges = []
        self.neighbors = []

    def add(self, substrate):
        self.edges.append(substrate)

    def addNeighbor(self, node):
        pass
            
        

