import math, sys, random

import pygame
from pygame.locals import *
from pygame.color import *

import pymunk as pm
from pymunk import Vec2d
import pymunk.pygame_util

class cell(object):
    """
    Every attribute of a cell will be written here.
    """
    def __init__(self,space,num,pos,radius=3,mass=0.1):
        """
        pos should be a tuple (x,y)
        num is the cell ID
        No two cells should share the same num.
        These are initial attributes. For updated attributes, use call functions
        """
        ## Cell Attributes
        position=pos
        self.num=num #Cell id
        radius=3 #Radius of hemisphere at cell end
        length=radius*2 #The length of the cell measure from the center of the two hemisphere
        height=radius*2 #The height of the box connecting the two hemisphere. = diameter
        mass=mass
        ## Initializing cell in space
        self.space=space
        inertia=pymunk.moment_for_poly(mass,\
        [(0,0), (0,height), (length,height), (height,0)], (0,0))
        self.body=pm.Body(mass,inertia)
        body.position = pos
        box = pm.Poly(body, [(0,0), (0,height), (length,height), (height,0)],\
        (-length/2.,-height/2.)) # Connecting box
        left = pm.Circle(body,radius,(-length/2.,0)) # left hemisphere
        right = pm.Circle(body,radius,(length/2.,0)) # right hemisphere
        left.color = pg.color.THECOLORS["red"]
        right.color = pg.color.THECOLORS["red"]
        box.color = pg.color.THECOLORS["red"]
        space.add(self.body, box, left, right)
        
    def __getattr__(self, attr):
        """
        Gets basic attribute
        """
    def get_pos():
        """
        Gets the position of the cell
        """
        pass
    def get_neighbour():
        """
        Depends on distance threshold
        Returns a list of tuple (neighbouring cell id, distance from neighbour)
        """
        pass
    def get_parents():
        """
        Gets a list of all parent cell id
        """
        pass
        
class cellspace():
    """
    Keeps track of all cells. Decides which cells are active, which are not active. 
    """
    def __init__(self):
        space = pm.Space()
        self.cells=[]
        self.active_cells=[]
        self.

        
    def add_cell():
        self.cells=()
    

def add_cell():


