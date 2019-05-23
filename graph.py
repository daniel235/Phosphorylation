import pandas as pd


class Graph:
    '''
        This class creates the graph using the weight of edges
        as # of substrates kinases share.

        graph class includes nodes edges

        node class uses name it's edges and neighbors

        edge class uses substrates names and number of edges.
    '''
    def __init__(self):
        self.graph = []
        self.kinaseNodes = []
        self.kinaseNodesNames = []
        self.kinaseData = pd.read_csv("./data/Kinase_Substrates.txt", delimiter="\t")
        self.edges = []
        self.edgeNames = []


    def createNode(self, name):
        node = Node(name)
        self.kinaseNodes.append(node)
        self.kinaseNodesNames.append(node.name)
        return node


    def add(self, node):
        self.graph.append(node)


    def addEdge(self, substrate):
        self.edgeNames.append(substrate)

        #create edge
        e = Edge(substrate)
        self.edges.append(e)

        return e 

    
    def grabEdge(self, substrate):
        if substrate in self.edgeNames:
            index = self.edgeNames.index(substrate)
            edge = self.edges[index]

        else:
            edge = self.addEdge(substrate)

        return edge



    def createGraph(self):
        for i in range(len(self.kinaseData)):
            #grab substrate object
            s = self.kinaseData[i][1]
            s = self.grabEdge(s)

            #check if kinase in nodes
            if self.kinaseData["Kinase"][i] in self.kinaseNodesNames:
                nodeIndex = self.kinaseNodesNames.index(self.kinaseData["Kinase"][i])
                node = self.kinaseNodes[nodeIndex]

            else:
                node = self.createNode(self.kinaseData["Kinase"][i])

            #grab nodes associated with this substrate
            subNodes = s.nodes

            #create neighbors or add neighbors to this node
            for i in range(len(subNodes)):
                if subNodes[i] in node.neighbors:
                    node.neighbors[subNodes[i]] += 1
                else:
                    node.neighbors[subNodes[i]] = 1

                
        



class Node:
    def __init__(self, name):
        self.name = name
        self.edges = []
        self.neighbors = {}

    def add(self, substrate):
        edge = Edge(substrate)
        self.edges.append(edge)

            
        

class Edge:
    def __init__(self, name):
        self.name = name
        self.nodes = []

    #add node object 
    def addNode(self, node):
        self.nodes.append(node)
