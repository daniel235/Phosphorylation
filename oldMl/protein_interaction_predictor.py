import numpy as np
import pandas as pd
import pygame


class protein_node():
    def __init__(self, name):
        self.name = name
        self.protein_edges = []


class protein_interaction_net():
    def __init__(self, objects):
        self.protein_names = []
        self.protein_objects = objects
        self.protein_interaction_data = np.array(pd.read_csv('./data/Protein_Protein_Interaction.txt', delimiter="\t"))
        self.nodes = []

    #might create gui for protein network
    def create_gui(self):
        pass

    #add edge
    def add_edge(self, node_object, edge):
        #edge has to go both ways (no duplicates)!
        node_object.protein_edges.append(edge)

    #create network
    def network(self):
        for i in range(len(self.protein_interaction_data)):
            first = self.protein_interaction_data[i][0]
            second = self.protein_interaction_data[i][1]

            #grab node objects from node array
            first_node = self.grab_node(first)
            second_node = self.grab_node(second)

            #check for each node in edge arrays
            if first_node not in second_node.protein_edges and second_node not in first_node.protein_edges:
                first_node.protein_edges.append(second_node)

        print(self.nodes[0].protein_edges[0])



    def create_node(self, name):
        n = protein_node(name)
        self.nodes.append(n)
        return n


    def grab_node(self, name):
        '''for i in range(len(self.nodes)):
            if self.nodes[i] == name:
                return self.nodes[i]'''

        try:
            return self.nodes.index(name)

        except:
            return self.create_node(name)

