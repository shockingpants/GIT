import math, sys, random

import pygame as pg
import pygame.locals as pl
import pygame.color as pc

import pymunk as pm
from pymunk import Vec2d
import pymunk.pygame_util

#########################################
#####		Initializing
#########################################
pg.init()
screen = pg.display.set_mode((600, 600))
clock = pg.time.Clock()
running = 1
balls = []
space = pm.Space()
#space.gravity = (0.0, -900.0)
space.add_collision_handler(0, 0, None, None, None, None, surface=screen)
time=0
dt = 1.0/60.0
#########################################
#####		Helping function
#########################################
def to_pygame(p):
	"""Small hack to convert pymunk to pygame coordinates"""
	return int(p.x), int(-p.y+600)

#########################################
#####		    Run
#########################################
#radius=25
#mass=0.1
#inertia = pm.moment_for_circle(mass,0, radius, (0,0))
#body = pm.Body(mass, inertia)
#body.position = 200,400
#shape = pm.Circle(body,radius,(0,0))
#shape.friction = 0.5
#shape.elasticity = 0.5
#body2 = pm.Body(mass, inertia)
#body2.position = 100,300
#shape2 = pm.Circle(body2,radius,(0,0))
#shape2.friction = 0.5
#shape2.elasticity = 0.5
boxw=60
boxh=30
mass=0.1
inertia=pymunk.moment_for_poly(0.1,[(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (0,0))
#inertia=pymunk.moment_for_poly(0.1,[(0,0), (0,1), (1,1), (1,0)], (-boxw/2.,-boxh/2.))
body = pm.Body(mass, inertia)
body.position = 200,150
shape1 = pm.Circle(body,boxh/2.,(-boxw/2.,0))
shape2 = pm.Circle(body,boxh/2.,(boxw/2.,0))
shape3 = pm.Poly(body, [(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (-boxw/2.,-boxh/2.))
shape1.color = pg.color.THECOLORS["red"]
shape2.color = pg.color.THECOLORS["red"]
shape3.color = pg.color.THECOLORS["red"]
#shape3 = pm.Poly(body, [(0,0), (0,boxh), (boxw,boxh), (boxw,0)],(boxw/2.,boxh/2.))
space.add(body, shape3, shape2, shape1)

bodya = pm.Body(mass, inertia)
bodya.position = 100,200
bodya.angle = math.pi/2
shape1a = pm.Circle(bodya,boxh/2.,(-boxw/2.,0))
shape2a = pm.Circle(bodya,boxh/2.,(boxw/2.,0))
shape3a = pm.Poly(bodya, [(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (-boxw/2.,-boxh/2.))
shape1a.color = pg.color.THECOLORS["red"]
shape2a.color = pg.color.THECOLORS["red"]
shape3a.color = pg.color.THECOLORS["red"]
#shape3 = pm.Poly(body, [(0,0), (0,boxh), (boxw,boxh), (boxw,0)],(boxw/2.,boxh/2.))
space.add(bodya, shape3a,shape2a,shape1a)



while running:
    ##Events
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
        elif event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE:
            running = False
        elif event.type == pg.KEYDOWN and event.key == pg.K_p:
            pg.image.save(screen, "contact_with_friction.png")
    #radius+=0.2
    #shape.unsafe_set_radius(radius)
    boxw+=1

    inertia=pymunk.moment_for_poly(0.1,[(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (0,0))
    body.moment=inertia
    bodya.moment=inertia

    shape3.unsafe_set_vertices([(0,0), (0,boxh), (boxw,boxh), (boxw,0)],(-boxw/2.,-boxh/2.))
    shape1.unsafe_set_offset((-boxw/2.,0))
    shape2.unsafe_set_offset((boxw/2.,0))

    shape3a.unsafe_set_vertices([(0,0), (0,boxh), (boxw,boxh), (boxw,0)],(-boxw/2.,-boxh/2.))
    shape1a.unsafe_set_offset((-boxw/2.,0))
    shape2a.unsafe_set_offset((boxw/2.,0))

    ## Draw
    screen.fill(pg.color.THECOLORS["black"])
    pm.pygame_util.draw(screen, space)
    
    ### Update physics
    space.step(dt)

    clock.tick(50)
    pg.display.set_caption("fps: " + str(clock.get_fps()))
    pg.display.flip()
    ## Misc
    time+=dt

#if __name__ == '__main__':
#    sys.exit()
