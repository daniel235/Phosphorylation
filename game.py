import pygame
import os, sys
import pickle


class gameManager():
    def __init__(self):
        #import data from pickles
        #get cluster groups
        with open("./data/pickles/clusterNodes", 'rb+') as f:
            self.clusterGroups = pickle.load(f)
        board = Surface()
        self.clusters = []

    def draw_objects(self):
        count = 0
        for cluster in self.clusterGroups[1]:
            self.clusters.append(Cluster(cluster.name))
            self.clusters[count].edges = cluster.edges
            self.clusters[count].data = cluster.data



    def update(self):
        while 1:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    sys.exit()

            

class Surface():
    def __init__(self):
        pygame.init()
        pygame.display.set_mode([800, 800])
        self.screen = pygame.display.get_surface()

    def draw(self, obj):
        pygame.draw(self.screen, obj.color, pygame.math.Vector2(obj.x, obj.y), obj.width)


class Cluster():
    def __init__(self, name):
        self.edges = {}
        self.name = name
        self.width = 20
        self.color = None
        self.x = None
        self.y = None
        self.data = None


gm = gameManager()
gm.update()