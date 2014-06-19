import sys
import pygame
from pygame.locals import *
from pygame.color import *
import math

############################
width = height = 600
### PyGame init
pygame.init()
screen = pygame.display.set_mode((width,height))
clock = pygame.time.Clock()
running = True
font = pygame.font.SysFont("Arial", 16)
### Clear screen
screen.fill(pygame.color.THECOLORS["black"])
### Draw stuff
pygame.draw.arc(screen, (255,0,0,255), [(0,0),(0,1),(1,1),(1,0)], math.pi/2, math.pi*3/2, 4)
pygame.display.flip()
fps=50
clock.tick(fps)
