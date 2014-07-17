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

from types import MethodType
instance2 = A_Class()
setattr(instance, fn.__name__, MethodType(fn, instance, type(instance)))

class biochemistry(object):
	"""
	This uses the SDE integrator. 
	Species are in the form of concentration
	May eventually include sbml, or may create a separate class for it
	Will want to add symbolic py eventually.
	Assumes two compartments. Cell and solution
	"""
	def __init__(self,func,y0,t0=0,dt=0.1,vol=None,param=None,integrator="SDE",method="rk"):
		##{{{
		"""
		func is the function coding the ODE
		vol is volume. It is also a proxy to vary the noise
		param is the parameter for func
		integrator is the type of integrator. SDE or ODE
		method is either "euler" (Euler M) or "rk" (Runge Kutta)
		"""
		self.t1=t1
		self.t0=t0
		self.time=[t0]
		self.dt=dt

		self.y0=np.array(y0)
		#Should insert an assert statement that checks func and y0
		self.values=np.array([y0])

		self.cell=cell #Cell
		self.cs=cell.cellspace #Cellspace
		self.sols=cell.cellspace.solspace #Solspace

		if integrator=="SDE":
			self.de=SDE_integrator(func,method=method)
			self.de.set_initial_value(y0,t0) #Adds parameter
			self.de.set_parameters(param)

		elif integrator=="ODE":
			self.de=od(func).set_integrator(integrator,method=int_method,**kwargs) 
			self.de.set_initial_value(y0,t0).set_f_params(param) #Adds parameter
		##}}}

	def assign_cell(self,cell):
		##{{{
		"""
		Assigns a cell to the biochemistry class.
		Did not add this to __init__ cuz user will want to initialize this before the cell is created
		"""
		self.cell=cell
		##}}}

	def run_de(self,dt):
		##{{{
		"""
		helper function for running differential equation by a fixed timestep
		"""
		# Run DE
		if de.successful():
			de.integrate(de.t+dt)
			self.time.append(de.t)
			self.values=np.concatenate((values,np.array([de.y])),axis=0)
		else:
			print de.successful()
		##}}}
	
	def updatecells(self,dt):
		##{{{
		"""
		User defined method. Bind this method to the biochemistry instance
		User will have to define the method in the run script
		http://dietbuddha.blogspot.com/2012/12/python-metaprogramming-dynamically.html
		Should be used to change system properties based on biochemical state of cell
		"""
		pass
		##}}}
	
	def transport(self,dt):
		##{{{
		"""
		Method used to exchange molecules/concentration between cell and media
		"""
		pass
		##}}}

	def reporter(self,fmin,fmax,color=):
		##{{{
		"""
		Change cell color based on protein concentration
		"""
		self.cell.color=
		pass
		##}}}

	def run(self,dt):
		##{{{
		"""
		"""
		self.run_ode(dt)
		self.transport(dt)
		self.reporter(dt)
		self.updatecells(dt)
		##}}}

	def divide(self,cell):
		##{{{
		"""
		https://docs.python.org/2/library/copy.html
		"""
		pass
		##}}}
		
