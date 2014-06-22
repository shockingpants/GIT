import cellsim as cs
import math
import pygame as pg

space=cs.cellspace(dt=0.1)
space.add_cell(cs.cellp,(300,300),radius=4,mass=0.01)
#space.add_cell((200,240),angle=math.pi/2)
#pg.draw.polygon(space.screen, (255,0,0,0), C.box.get_vertices(), 0)
#pg.display.flip()
space.run()
#space.stop()
