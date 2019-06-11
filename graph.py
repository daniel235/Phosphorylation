import pandas as pd
import cluster_data
#create graph data for cytoscape
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
        e = SubEdge(substrate)
        self.edges.append(e)

        return e 

    
    def grabEdge(self, substrate):
        if substrate in self.edgeNames:
            index = self.edgeNames.index(substrate)
            edge = self.edges[index]

        else:
            edge = self.addEdge(substrate)

        return edge


    def connectGraph(self):
        nodesStillNeeded = self.kinaseNodesNames

        #create sif file
        filename = "./graphMod.sif"
        with open(filename, 'a') as f:
            for i in range(len(nodesStillNeeded)):
                if not nodesStillNeeded:
                    return 
                #create edge between nodes neighbors 
                #grab node
                index = self.kinaseNodesNames.index(nodesStillNeeded[i])
                node = self.kinaseNodes[index]
                newLen = len(nodesStillNeeded)

                print("before loop")
                print(node.name)
                print(len(node.neighbors))
                for neigh, value in node.neighbors.items():
                    print(neigh)
                    #add edges 
                    e = Edge(node, neigh)
                    e.weight = value 

                    #write to file
                    if(e.weight != 1):
                        f.write(e.pair[0].name + " " + e.pair[1].name + " " + str(e.weight) + "\n")

                    #remove from to do list / assuming neighbors are node objects(not sure)
                    #nodesStillNeeded.remove(neigh.name)

        #!todo cutoff
        #plot histogram of kinase to cutoff
    
    def createGraph(self):
        for i in range(len(self.kinaseData)):
            #grab substrate object
            s = self.kinaseData["Substrate"][i]
            print("substrate " + s)
            s = self.grabEdge(s)
            print(s.nodes)


            #check if kinase in nodes
            if self.kinaseData["Kinase"][i] in self.kinaseNodesNames:
                nodeIndex = self.kinaseNodesNames.index(self.kinaseData["Kinase"][i])
                node = self.kinaseNodes[nodeIndex]

            else:
                node = self.createNode(self.kinaseData["Kinase"][i])


            #add node to this sub edge
            s.nodes.append(node)
            #grab nodes associated with this substrate
            subNodes = s.nodes

            #create neighbors or add neighbors to this node
            for i in range(len(subNodes)):
                if subNodes[i] in node.neighbors and subNodes[i] != node:
                    node.neighbors[subNodes[i]] += 1
                else:
                    node.neighbors[subNodes[i]] = 1

        print("connect graph now")
        self.connectGraph()


    def createKinaseClassHistograms(self, n):
        c = cluster_data.ClusterData()
        c.set_arbitrary_kinase_class(n)
        wKinase = c.weakKinase
        rKinase = c.strongKinase
        
        return wKinase, rKinase



class Node:
    def __init__(self, name):
        self.name = name
        self.edges = []
        self.neighbors = {}

    def add(self, substrate):
        edge = SubEdge(substrate)
        self.edges.append(edge)

            
        

class SubEdge:
    def __init__(self, name):
        self.name = name
        self.nodes = []

    #add node object 
    def addNode(self, node):
        self.nodes.append(node)


class Edge:
    def __init__(self, node1, node2):
        self.pair = (node1, node2)
        self.weight = 0

