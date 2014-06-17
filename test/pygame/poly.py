import math, sys, random
import pygame as pg
import pygame.locals as pl
import pygame.color as pc
import pymunk as pm
from pymunk import Vec2d
import pymunk.pygame_util
#import pdb
#########################################
#####		Initializing
#########################################
pg.init()
screen = pg.display.set_mode((600, 600))
clock = pg.time.Clock()
running = 1
space = pm.Space()
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
boxw=60
boxh=30
mass=0.1
inertia=pymunk.moment_for_poly(0.1,[(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (0,0))

body = pm.Body(mass, inertia)
body.position = 200,100
shape3 = pm.Poly(body, [(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (-boxw/2.,-boxh/2.))
shape3.color = pg.color.THECOLORS["red"]
space.add(body, shape3)

bodya = pm.Body(mass, inertia)
bodya.position = 200,200
bodya.angle = math.pi/2+math.pi
shape3a = pm.Poly(bodya, [(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (-boxw/2.,-boxh/2.))
shape3a.color = pg.color.THECOLORS["red"]
space.add(bodya, shape3a)



while running:
    ##Events
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
        elif event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE:
            running = False
        elif event.type == pg.KEYDOWN and event.key == pg.K_p:
            pg.image.save(screen, "contact_with_friction.png")

    boxw+=0.2
    inertia=pymunk.moment_for_poly(0.1,[(0,0), (0,boxh), (boxw,boxh), (boxw,0)], (0,0))
    body.moment=inertia
    bodya.moment=inertia
    shape3.unsafe_set_vertices([(0,0), (0,boxh), (boxw,boxh), (boxw,0)],offset=(-boxw/2.,-boxh/2.))
    #shape3.offset=(-boxw/2.,-boxh/2.)
    shape3a.unsafe_set_vertices([(0,0), (0,boxh), (boxw,boxh), (boxw,0)],offset=(-boxw/2.,-boxh/2.))
    #shape3a.offset=(-boxw/2.,-boxh/2.)
    #raw_input()

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

