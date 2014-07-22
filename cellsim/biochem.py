###############################################################################################################
# CellSim
#
# Copyright Jonathan Teo   
# Contact   jonteo@mit.edu
# 2014
###############################################################################################################
import numpy as np
import math, sys, random
import copy

import pygame as pg
import pygame.locals as pl
import pygame.color as pc
import pymunk as pm

import fipy as fp
import fipy.tools.numerix as fnumerix

import scipy
import scipy.integrate
od=scipy.integrate.ode

from cellsim.integrator import *
from cellsim.viewer import hex_to_rgb

#from types import MethodType
#instance2 = A_Class()
#setattr(instance, fn.__name__, MethodType(fn, instance, type(instance)))

class biochemistry(object):
	"""
	This uses the SDE integrator. 
	Species are in the form of concentration
	May eventually include sbml, or may create a separate class for it
	Will want to add symbolic py eventually.
	Assumes two compartments. Cell and solution
	"""
	def __init__(self,func,names,y0,t0=0,dt=0.1,vol=None,param=None,integrator="SDE",method="rk"):
		##{{{
		"""
		func is the function coding the ODE
		vol is volume. It is also a proxy to vary the noise
		param is the parameter for func
		integrator is the type of integrator. "SDE" or "ODE"
		method is either "euler" (Euler M) or "rk" (Runge Kutta)
		"""
		self.t0=t0
		self.time=[t0]
		self.dt=dt
		self.y0=np.array(y0)
		self.integrator=integrator
		self.method=method
		self.param=param
		self.func=func
		self.vol=vol
		assert len(names)==len(y0)
		self.names=names #These are the names of the species. The order should correspond to order of
						 # appearance in func.
		#Should insert an assert statement that checks func and y0
		self.values=np.array([y0])
		self.de=self.set_integrator(self.t0,self.y0,self.integrator)
		##}}}

	def set_integrator(self,t,y,integrator="SDE"):
		##{{{
		"""
		"""
		if integrator=="SDE":
			de=SDE_integrator(self.func,vol=self.vol,method=self.method)
			de.set_initial_value(y,t) #Adds parameter
			de.set_parameters(self.param)

		elif integrator=="ODE":
			ODEintegrator="vode"
			de=od(self.func).set_integrator(self.ODEintegrator,method=method,**kwargs) 
			de.set_initial_value(y,t).set_f_params(self.param) #Adds parameter

		return de

		##}}}

	def assign_cell(self,cell):
		##{{{
		"""
		Assigns a cell to the biochemistry class.
		Did not add this to __init__ cuz user will want to initialize this before the cell is created
		"""
		self.cell=cell
		self.cs=cell.cellspace #Cellspace
		self.sols=cell.cellspace.solspace #Solspace
		##}}}

	def run_de(self,dt):
		##{{{
		"""
		helper function for running differential equation by a fixed timestep
		"""
		# Run DE
		de=self.de
		if de.successful():
			de.integrate(de.t+dt)
			self.time.append(de.t)
			self.values=np.concatenate((self.values,np.array([de.y])),axis=0)
			#self.time=de.t
			#self.values=de.y
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

	def reporter(self,ind,color="#FF0000",fmin=0,fmax=1):
		##{{{
		"""
		Change cell color based on protein concentration
		ind indicates which variable
		vvvvaaaaaaaaa
		fmin and fmax will be changed to 
		"""
		if type(color) is str:
			color = color.lstrip('#')
			lv = len(color)
			assert lv==6
			RGB=tuple(int(color[i:i+lv/3], 16) for i in range(0, lv, lv/3))
		elif type(color) is tuple:
			assert len(color)==3
			RGB=color
		else:
			raise TypeError('color must be hex string or RGB tuple.')
		self.reportercolor=RGB
		self.ind=ind
		self.fmin=fmin
		self.fmax=fmax
		##}}}

	def _update_color(self):
		##{{{
		"""
		Update color of cell based on reporter
		"""
		if "reportercolor" in vars(self):
			ratio=(self.values[-1,self.ind]-self.fmin)/(self.fmax-self.fmin)
			if ratio < 0.0:
				self.cell.color=(255,255,255)
			elif ratio > 1.0:
				self.cell.color=self.reportercolor
			else:
				self.cell.color=np.array((255,255,255))-ratio*(np.array((255,255,255))-np.array(self.reportercolor))
		else:
			pass

		##}}}
	
	def get_latest(self,name):
		##{{{
		"""
		Get latest value of species corresponding to the name
		"""
		return self.values[self.names.index[name]][-1]
		##}}}

	def run(self,dt):
		##{{{
		"""
		"""
		self.run_de(dt)
		self.transport(dt)
		self._update_color()
		self.updatecells(dt)
		##}}}
	
	def divide(self):
		##{{{
		"""
		Need to determine which variable in biochem should be copied over
		Use copy.copy
		https://docs.python.org/2/library/copy.html
		"""
		cp=copy.copy #To cope instead of referencing a variable
		temp=biochemistry(cp(self.func),cp(self.names),cp(self.de.y),t0=cp(self.de.t),dt=cp(self.dt),vol=cp(self.vol),param=cp(self.param),\
		integrator=cp(self.integrator),method=cp(self.method))
		temp.de=temp.set_integrator(cp(self.de.t), cp(self.de.y), cp(self.integrator))
		temp.values=cp(self.values)
		temp.time=cp(self.time)
		temp.reportercolor=cp(self.reportercolor)
		temp.ind=cp(self.ind)
		temp.fmin=cp(self.fmin)
		temp.fmax=cp(self.fmax)
		return temp

		##}}}
	
		##}}}
