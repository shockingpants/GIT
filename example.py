import cellsim as cs
import math
import pygame as pg

space=cs.cellspace(dt=0.1)
b=space.add_cell(cs.cellp,(300,300),radius=4,mass=0.01)
c=space.add_cell(cs.cellp,(200,300),radius=4,mass=0.01)
space.solspace.add_species('sig1',degradation=1, diffusion=1)
#space.solspace.add_species('sig2',degradation=1, diffusion=1,color="#FF0000")
space.run()



#space.add_cell((200,240),angle=math.pi/2)
#pg.draw.polygon(space.screen, (255,0,0,0), C.box.get_vertices(), 0)
#pg.display.flip()
#space.stop()
