import pygame
import os, sys
import pickle


class gameManager():
    def __init__(self):
        #import data from pickles
        #get cluster groups
        with open("./data/pickles/clusterNodes", 'rb+') as f:
            self.clusterGroups = pickle.load(f)
        self.board = Surface()
        self.clusters = []

    def draw_objects(self):
        pos = (50, 50)
        count = 0
        for cluster in self.clusterGroups[1]:
            self.clusters.append(Cluster(cluster.name))
            self.clusters[count].edges = cluster.edges
            self.clusters[count].data = cluster.data
            #set position for cluster width(20)
            self.clusters[count].x = int((count + 1) * pos[0])
            self.clusters[count].y = int(pos[1])
            self.clusters[count].pos = pygame.math.Vector2(self.clusters[count].x, self.clusters[count].y)
            count += 1
        
        #after finish initializing cluster properties draw clusters
        for cluster in self.clusters:
            self.board.draw(cluster)


    def update(self):
        while 1:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    sys.exit()
                
                else:
                    pygame.display.update()

            

class Surface():
    def __init__(self):
        pygame.init()
        pygame.display.set_mode([800, 800])
        self.screen = pygame.display.get_surface()

    def draw(self, obj):
        print(type(obj.color), type(obj.pos), type(obj.width))
        pygame.draw.circle(self.screen, obj.color, (int(obj.pos.x), int(obj.pos.y)), int(obj.width))


class Cluster():
    def __init__(self, name):
        self.edges = {}
        self.name = name
        self.width = 20
        self.color = pygame.Color(20, 100, 100)
        self.x = None
        self.y = None
        self.pos = None
        self.data = None


gm = gameManager()
gm.draw_objects()
gm.update()