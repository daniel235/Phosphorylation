import pygame
import compare_clusters
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
        self.font = pygame.font.get_default_font()


    def draw_objects(self):
        pos = (60, 60)
        count = 0
        countx = 0
        black = (0, 0, 0)
        for subtype in range(len(self.clusterGroups)-1):
            for cluster in self.clusterGroups[subtype]:
                #write text
                #text = self.font.render(cluster.name, True, black)
                #TextSurf = text
                #TextRect =  text.get_rect()
                
                
                self.clusters.append(Cluster(cluster.name))
                self.clusters[count].edges = cluster.edges
                self.clusters[count].data = cluster.data
                color = []
                for i in range(3):
                    self.clusters[count].color[i] = self.clusters[count].color[i] + (subtype * 50)
        
                self.clusters[count].color = pygame.Color(self.clusters[count].color[0], self.clusters[count].color[1], self.clusters[count].color[2])
                #set position for cluster width(20)
                self.clusters[count].x = int((countx + 1) * pos[0])
                self.clusters[count].y = int(pos[1]) + (subtype * 200)
                #TextRect.center = ((self.clusters[count].x),(self.clusters[count].y))
                #self.board.disp.blit(TextSurf, TextRect)
                self.clusters[count].pos = pygame.math.Vector2(self.clusters[count].x, self.clusters[count].y)
                count += 1
                countx += 1
        
            countx = 0
            
        #after finish initializing cluster properties draw clusters
        for cluster in self.clusters:
            self.board.draw(cluster)

        
        cnt = 0
        for cluster in self.clusters:
            cnt += 1
            #create all edged in each cluster
            for edge, score in cluster.edges.items():
                print(edge)
                #find cluster 
                for famcluster in self.clusters:
                    if famcluster.name == edge:
                        if score > .05:
                            color = (255, 255, 255)
                        else:
                            if score != 0:
                                color = ((200 - cnt), (cnt * 10), (cnt * 3))
                                self.board.drawLine(cluster, famcluster, color)
                        


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
        self.disp = pygame.display.set_mode([800, 800])
        self.screen = pygame.display.get_surface()


    def draw(self, obj):
        print(type(obj.color), type(obj.pos), type(obj.width))
        pygame.draw.circle(self.screen, obj.color, (int(obj.pos.x), int(obj.pos.y)), int(obj.width))


    def drawLine(self, obj, obj2, color):
        pygame.draw.lines(self.screen, color, True, [(obj.x, obj.y), (obj2.x, obj2.y)])



class Cluster():
    def __init__(self, name):
        self.edges = {}
        self.name = name
        self.width = 30
        self.color = [20, 100, 100]
        self.x = None
        self.y = None
        self.pos = None
        self.data = None

    

gm = gameManager()
gm.draw_objects()
gm.update()