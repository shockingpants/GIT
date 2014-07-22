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
space=cs.cellspace(dt=0.05)
names=['c','i']
bc1=cs.biochemistry(test_osc,names,y0=[2,2],vol=1e-22,param=param, integrator="SDE",method="rk")
#bc1.reporter('x',color=) Implement symbolic py for this
bc1.reporter(1,color="#00ff00")
b=space.add_cell(cs.cellp,(300,300),radius=4,length=15,mass=0.01,biochem=bc1)
#bc1=cs.biochem(test_osc,y0,vol=1e-21,param=param, integrator="SDE",method="rk")
#c=space.add_cell(cs.cellp,(200,300),radius=4,mass=0.01,biochem=bc1)
space.solspace.add_species('c',degradation=1, diffusion=10, kin=1, kout=1,datamax=0.1)
space.solspace.add_species('i',degradation=1, diffusion=1, kin=1, kout=1, color="#FF0000",datamax=0.1)
space.run()



#space.add_cell((200,240),angle=math.pi/2)
#pg.draw.polygon(space.screen, (255,0,0,0), C.box.get_vertices(), 0)
#pg.display.flip()
#space.stop()
