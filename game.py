import pygame
import os, sys



class gameManager():
    def __init__(self):
        #import data from pickles
        pass


    def draw_objects(self):
        pass

class Surface():
    def __init__(self):
        screen = pygame.display.get_surface()
        self.area = screen.get_rect()
        self.rect.topleft = 10, 10

    def draw(obj):
        pygame.draw(Surface.screen, obj.color, pygame.math.Vector2(obj.x, obj.y)  self.width)


class cluster():
    def __init__(self, name):
        self.edges = {}
        self.name = name
        self.width = 20
        self.color = None
        self.x = None
        self.y = None

