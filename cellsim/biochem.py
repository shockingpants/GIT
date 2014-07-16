###############################################################################################################
# CellSim
#
# Copyright Jonathan Teo   
# Contact   jonteo@mit.edu
# 2014
###############################################################################################################
import numpy as np
import math, sys, random

import pygame as pg
import pygame.locals as pl
import pygame.color as pc
import pymunk as pm

import fipy as fp
import fipy.tools.numerix as fnumerix

import scipy
import scipy.integrate
od=scipy.integrate.ode

class biochemistry(object):
	"""
	This uses the SDE integrator. 
	Species are in the form of concentration
	May eventually include sbml, or may create a separate class for it
	Will want to add symbolicpy eventually
	"""
	def __init__(self,cell,func,y0,t0=0,t1=50,dt=0.1,param=None,integrator='SDE'):
		"""
		"""
		self.t1=t1
		self.t0=t0
		self.time=[t0]
		self.dt=dt

		self.y0=np.array(y0)
		#Should insert an assert statement that checks func and y0
		self.values=np.array([y0])

		self.cell=cell
		if integrator=="SDE":
			self.de=SDE_integrator(func,method=method)
			self.de.set_initial_value(y0,t0) #Adds parameter
			self.de.set_parameters(param)

		elif integrator=="ODE":
			self.de=od(func).set_integrator(integrator,method=int_method,**kwargs) 
			self.de.set_initial_value(y0,t0).set_f_params(param) #Adds parameter
		##}}}


	def run(self,dt):
		"""
		"""
		# Run ODE
		if de.successful():
			de.integrate(de.t+dt)
			self.time.append(de.t)
			self.values=np.concatenate((values,np.array([de.y])),axis=0)
		else:
			print de.successful()


	def divide(self,cellspace):
		"""
		https://docs.python.org/2/library/copy.html
		"""
		
