import cellsim as cs
import numpy as np
#import math
#import pygame as pg
#from types import MethodType
#instance = A_Class()
#setattr(instance, fn.__name__, MethodType(fn, instance, type(instance)))

def test_osc(t,y,param=None):
	##{{{
	c,i=y
	
	p=param
	beta=p['beta']
	kappa=p['kappa']

	dcdt=c**2/(beta*i**2+c**2)-c
	didt=kappa*(c-i)
	return np.array([dcdt,didt])
	##}}}
param=dict([("beta",4),("kappa",0.5)])
space=cs.cellspace(dt=0.1)
#bc1=cs.biochemistry(test_osc,y0=[2,2],vol=1e-21,param=param, integrator="SDE",method="rk")
#cs.reporter('x',color=) Implement symbolic py for this
#bc1.reporter(1,color="#ff0000")
b=space.add_cell(cs.cellp,(300,300),radius=4,mass=0.01)
#bc1=cs.biochem(test_osc,y0,vol=1e-21,param=param, integrator="SDE",method="rk")
#c=space.add_cell(cs.cellp,(200,300),radius=4,mass=0.01,biochem=bc1)
space.solspace.add_species('sig1',degradation=1, diffusion=10)
space.solspace.add_species('sig2',degradation=1, diffusion=1,color="#FF0000")
space.run()



#space.add_cell((200,240),angle=math.pi/2)
#pg.draw.polygon(space.screen, (255,0,0,0), C.box.get_vertices(), 0)
#pg.display.flip()
#space.stop()
