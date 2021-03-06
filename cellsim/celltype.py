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

class cell(object):
	##{{{
	"""
	Every attribute of a cell will be written here.
	Should not be called directly since it contains no shape information. Should call using subclasses
	"""
	def __init__(self,cellspace,num,pos,angle,biochem=None,**kwargs):
		##{{{
		"""
		pos should be a tuple (x,y)
		num is the cell ID
		No two cells should share the same num.
		These are initial attributes.
		"""
		## Cell Attributes
		self.cellspace=cellspace
		self.space=cellspace.space #pymunk space class
		self.num=num #Cell id
		## Initialize states
		self.biochem=biochem # this should be the biochem instance. 
		#NOTE: MAKE SURE EACH CELL REFERENCES A UNIQUE biochem instance. Creat a new instance for each cell!!!
		if biochem is not None:
			biochem.assign_cell(self) # Assigns cell to biochem instance 
		## Intialize other attributes
		for i in kwargs.iterkeys(): 
			vars(self)[i]=kwargs[i]
		##}}}
	def _update_states(self,dt):
		##{{{
		"""
		Upadtes self.states. Stores them
		"""
		if self.biochem is not None:
			self.biochem.run(dt)
		##}}}
	def _check_states(self):
		##{{{
		"""
		Depending on what the chemical state is, do sthg
		Most of this will be coded from the biochem class instead.
		In fact, it should be user defined. 
		"""
		# Change to grow, quienescence, death
		# Change from active list to inactive list
		pass
		##}}}
	def help(self):
		##{{{	
		"""
		To call for attributes from body, just use [self].body.position. Replace self with variable name
		"""
		pass
		##}}}
	def get_pos(self):
		##{{{
		"""
		Gets the position of the cell
		"""
		return self.body.position
		##}}}
	def get_neighbour(self):
		##{{{
		"""
		Depends on distance threshold
		Returns a list of tuple (neighbouring cell id, distance from neighbour)
		"""
		pass
		##}}}
	def get_parents(self):
		##{{{
		"""
		Gets a list of all parent cell id
		"""
		pass
		##}}}
	
	##}}}
class cellp(cell):
	##{{{
	"""
	Cell approximated by a polygon
	1 px is 0.25um	
	"""
	def __init__(self,cellspace,num,pos,angle,radius=5,length=None,vertices=6,mass=0.1,biochem=None,**kwargs):
		##{{{
		## Cell Attributes
		super(cellp, self).__init__(cellspace,num,pos,angle,biochem=biochem,**kwargs)
		self.radius=radius #Radius of hemisphere at cell end. Radius is in uM
		self.length=radius*2 if length is None else length #The length of the cell measured from the center of the two hemisphere
		self.height=self.radius*2
		self.vertices=vertices #Number of vertices in the hemisphere
		self.mass=mass #Self Explanatory
		self.color=pg.color.THECOLORS["red"]

		self.growthrate=0.05 #Growth rate can be a function of a state variable eventually. um/s growth in length
		self.divmean=4.*self.radius #Threshold for division
		self.divstd=0.1*self.divmean
		self.divthreshold=self.divmean+self.divstd*np.random.randn(1) #Threshold decided at the start to reduce computational cost

		self.cycle='grow'
		self.kwargs=kwargs	# This is done so that when the cell divides, any special attribute
							# will be transfered to the daughter cell.
		for i in kwargs.iterkeys(): 
			vars(self)[i]=kwargs[i]

		## Initializing cell in space
		self._genpoly()
		inertia=pm.moment_for_poly(self.mass,self.ver, (0,0))
		self.body=pm.Body(self.mass,inertia)
		self.body.position = pos
		self.body.angle = angle
		#-----------
		self.box = pm.Poly(self.body,self.ver,(-self.length/2.,-self.height/2.)) # Connecting box
		#-----------
		self.box.elasticity=0.5
		self.box.friction=0
		self.box.color = self.color
		self.space.add(self.body, self.box)
		##}}}
	def	_update_cell(self):
		##{{{
		"""
		Updates Inertia, Updates length, height etc for pygame
		Does not update position.
		"""
		self._genpoly()
		self.body.moment=pm.moment_for_poly(self.mass,self.ver, (0,0)) #Sets new inertia
		self.box.unsafe_set_vertices(self.ver,(-self.length/2.,-self.height/2.))
		##}}}
	def _check_division(self):
		##{{{
		if self.length>self.divthreshold:
			self.divide()
		##}}}
	def _genpoly(self):
		##{{{
		"""
		Generate Vertices of a polygon
		"""
		ang=math.pi/self.vertices
		rotmat=np.array([[math.cos(ang),-math.sin(ang)],[math.sin(ang),math.cos(ang)]])
		verright=[np.array([0,-self.radius])] #Vertices for right hemisphere
		verleft=[np.array([0,self.radius])]
		for i in xrange(self.vertices):
			verright.append(np.dot(rotmat,verright[-1]))
			verleft.append(np.dot(rotmat,verleft[-1]))
		#Translating hemisphere based on cell length
		verright=np.array(verright)
		verleft=np.array(verleft)
		verright=verright+[self.length,self.radius]
		verleft=verleft+[0,self.radius]
		self.ver=np.vstack(((0,0),verright,verleft[0:-1]))
		##}}}
	def draw(self):
		##{{{
		"""
		Draw with py.games
		"""
		screen=self.cellspace.viewer.screen
		##----- Box------
		ps = self.box.get_vertices()
		ps = [to_pygame(p, screen) for p in ps]
		ps += [ps[0]]
		pg.draw.polygon(screen, self.color, ps, 0) #Jon
		pg.draw.polygon(screen, (0,0,0,255), ps, 1)
		#pg.draw.line(screen, (0,0,0,255), ps[1],ps[2], 2) #Jon
		#pg.draw.line(screen, (0,0,0,255), ps[3],ps[0], 2) #Jon
		##}}}
	def evolve(self,dt):
		##{{{
		"""
		Helper function to decide what happens to the cell
		Move cell forward by one time step
		3 states. grow, quie, dead
		"""
		self._update_states(dt)
		if self.cycle=='grow':	
			self.grow(dt)
			self._check_division()
		elif self.cycle=='quie':
			pass
		elif self.cycle=='dead':
			pass
		##}}}
	def divide(self):
		##{{{
		"""
		Triggered when cell is primed to divide
		"""
		###### Divides states into two
		#statea,stateb=self.split_state()
		# This is for SDE. No Change in concentration
		if self.biochem is not None:
			bc=self.biochem.divide()
		else:
			bc=None
		###### Gets position of daughter cells
		oldpos=self.body.position
		a=self.box.get_vertices()[5] #top left hand corner
		b=self.box.get_vertices()[13]
		c=self.box.get_vertices()[12]
		d=self.box.get_vertices()[6]
		ad=d-a #Vector 
		ab=b-a

		ran=np.random.randn()*0.05**2
		pos1=(a+ad*0.5*(0.5-(self.radius/self.length))+(ab/2.0))#*(1+ran) #a+ad direction+ab direction
		pos2=(d-ad*0.5*(0.5-(self.radius/self.length))+(ab/2.0))#*(1-ran)
		#Update Primary Daughter
		self.length=(self.length/2.0-self.radius)
		self.body.position=pos1
		self._update_cell()
		#Updating New Daughter
		self.kwargs
		temp=self.cellspace.add_cell(cellp,pos2,self.body.angle,length=self.length,radius=self.radius,biochem=bc,**self.kwargs) #TBDadd state next time
		temp.body.velocity=self.body.velocity
		temp.body.angular_velocity=self.body.angular_velocity
		temp.color=self.color
		# TAKENOTE, MAKE SURE THERE ARE NO NEW KEYWORDS ADDED, OTHERWISE **KWARGS WILL INCREASE WITH EVERY DIVISION
		# AND WE DON'T WANT THAT	

		##}}}
	def grow(self,dt):
		##{{{
		"""
		Run to grow a cell
		"""
		self.length*=1+dt*self.growthrate #Assumes exponential growth for now. TBD: Check cell doubling time
		self._update_cell()
		##}}}

	##}}}
class cellc(cell):
	##{{{
	"""
	Cell approximated by a box and two circles
	"""
	def __init__(self,cellspace,num,pos,angle,radius=4,mass=0.1,**kwargs):
		##{{{
		## Cell Attributes
		super(cellc, self).__init__(cellspace,num,pos,angle,**kwargs)
		self.radius=radius #Radius of hemisphere at cell end
		self.length=self.radius*2 #The length of the cell measured from the center of the two hemisphere
		self.height=self.radius*2 #The height of the box connecting the two hemisphere. = diameter
		self.mass=mass #Self Explanatory
		self.color="red"
		self.growthrate=0.5 #Growth rate can be a function of a state variable eventually	
		self.cycle='grow'
		for i in kwargs.iterkeys(): 
			vars(self)[i]=kwargs[i]

		## Initializing cell in space
		inertia=pm.moment_for_poly(self.mass,\
		[(-self.length/2.0,-self.height/2.0), (-self.length/2.0,self.height/2.0), (self.length/2.0,self.height/2.0), (self.length/2.0,-self.height/2.0)],\
		(0,0))
		self.body=pm.Body(self.mass,inertia)
		self.body.position = pos
		self.body.angle = angle
		#-----------
		self.box = pm.Poly(self.body, [(0,0), (0,self.height), (self.length,self.height), (self.length,0)],\
		(-self.length/2.,-self.height/2.)) # Connecting box
		self.left = pm.Circle(self.body,self.radius,(-self.length/2.,0)) # left hemisphere
		self.right = pm.Circle(self.body,self.radius,(self.length/2.,0)) # right hemisphere
		#-----------
		self.left.color=self.right.color=self.box.color = pg.color.THECOLORS[self.color]
		self.space.add(self.body, self.left, self.right, self.box)
		##}}}
	def	_update_cell(self):
		##{{{
		"""
		Updates Inertia, Updates length, height etc.
		Does not update position.
		"""
		self.body.moment=pm.moment_for_poly(self.mass,\
		[(-self.length/2.0,-self.height/2.0), (-self.length/2.0,self.height/2.0), (self.length/2.0,self.height/2.0), (self.length/2.0,-self.height/2.0)],\
		(0,0)) #Sets new inertia
		self.box.unsafe_set_vertices([(0,0), (0,self.height), (self.length,self.height), (self.length,0)],\
		(-self.length/2.,-self.height/2.))
		self.left.unsafe_set_offset((-self.length/2.,0))
		self.right.unsafe_set_offset((self.length/2.,0))
		self.left.unsafe_set_radius(self.radius)
		self.right.unsafe_set_radius(self.radius)
		##}}}
	def _check_division(self):
		##{{{
		self.division=5*self.radius #Threshold for division
		if self.length>self.division:
			self.divide()
		##}}}
	def draw(self):
		##{{{
		"""
		Draw with py.games
		"""
		screen=self.cellspace.viewer.screen
		##-------Left-------
		circle_center = self.body.position + self.left.offset.rotated(self.body.angle)
		pl = to_pygame(circle_center, screen)
		color = pg.color.THECOLORS["red"]
		diam=int(self.radius)*2
		#pg.draw.arc(screen, (0,0,0,255), ps, math.pi/2, math.pi*3/2, 4)
		pg.draw.circle(screen, color, pl, int(self.radius), 0)
		pg.draw.circle(screen, color, pl, int(self.radius), 4)
		##-------Right-------
		circle_center = self.body.position + self.right.offset.rotated(self.body.angle)
		pr = to_pygame(circle_center, screen)
		diam=int(self.radius)*2
		#pg.draw.arc(screen, (0,0,0,255), ps, math.pi/2, math.pi*3/2, 4)
		pg.draw.circle(screen, color, pr, int(self.radius), 0)
		pg.draw.circle(screen, color, pr, int(self.radius), 4)
		##----- Box------
		ps = self.box.get_vertices()
		ps = [to_pygame(p, screen) for p in ps]
		ps += [ps[0]]
		pg.draw.polygon(screen, color, ps, 0) #Jon
		pg.draw.line(screen, (0,0,0,255), ps[1],ps[2], 2) #Jon
		pg.draw.line(screen, (0,0,0,255), ps[3],ps[0], 2) #Jon
		##}}}
	def evolve(self,dt):
		##{{{
		"""
		Helper function to decide what happens to the cell
		Move cell forward by one time step
		3 states. grow, quie, dead
		The biochecmical processed are always running. 
		"""
		self._update_states(dt)
		if self.cycle=='grow':
			self.grow(dt)
			self._check_division()
		elif self.cycle=='quie':
			pass
		elif self.cycle=='dead':
			pass
		##}}}
	def divide(self):
		##{{{
		"""
		Triggered when cell is primed to divide
		"""
		#statea,stateb=self.split_state()
		oldpos=self.body.position
		a,b,c,d=self.box.get_vertices()
		ad=d-a #Vector 
		ab=b-a
		ran=np.random.randn()*0.05**2
		pos1=(a+ad*0.5*(0.5-(self.radius/self.length))+(ab/2.0))*(1+ran) #a+ad direction+ab direction
		pos2=(d-ad*0.5*(0.5-(self.radius/self.length))+(ab/2.0))*(1-ran)
		#Update Primary Daughter
		self.length=(self.length/2.0-self.radius)
		self.body.position=pos1
		#self._update_cell()
		#Updating New Daughter
		temp=self.cellspace.add_cell(cellc,pos2,self.body.angle,length=self.length,radius=self.radius) #add state next time
		###### Divides states into two
		#statea,stateb=self.split_state()
		# This is for SDE. No Change in concentration
		#self.biochem.divide(temp)
		##}}}
	def grow(self,dt):
		##{{{
		"""
		Run to grow a cell
		"""
		self.length*=1+dt*self.growthrate #Assumes exponential growth for now. TBD: Check cell doubling time
		self._update_cell()
		#self._check_division()
		##}}}

	##}}}

#############################
####  Utility Function
#############################
def to_pygame(p, surface):
	##{{{
	"""Convenience method to convert pymunk coordinates to pygame surface 
	local coordinates
	"""
	flip_y=True
	if flip_y:
		return int(p[0]), surface.get_height()-int(p[1])
	else:
		return int(p[0]), int(p[1])
	##}}}
def from_pygame(p, surface):
	##{{{
	"""Convenience method to convert pygame surface local coordinates to 
	pymunk coordinates	  
	"""
	return to_pygame(p,surface)
	##}}}

