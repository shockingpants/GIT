import math, sys, random
import pygame as pg
import pygame.locals as pl
import pygame.color as pc
import pymunk as pm
from pymunk import Vec2d
import pymunk.pygame_util
import numpy as np
import pdb

class cell(object):
	##{{{
	"""
	Every attribute of a cell will be written here.
	"""
	def __init__(self,cellspace,num,pos,angle,radius=10,mass=0.1,**kwargs):
		##{{{
		"""
		pos should be a tuple (x,y)
		num is the cell ID
		No two cells should share the same num.
		These are initial attributes. For updated attributes, use call functions
		"""
		## Cell Attributes
		self.cellspace=cellspace
		self.num=num #Cell id
		self.radius=radius #Radius of hemisphere at cell end
		self.length=self.radius*2 #The length of the cell measured from the center of the two hemisphere
		self.height=self.radius*2 #The height of the box connecting the two hemisphere. = diameter
		self.mass=mass #Self Explanatory
		self.color="red"
		self.growthrate=0.5 #Growth rate can be a function of a state variable eventually
		self.division=10*self.radius #Threshold for division
		self.cycle='grow'
		## Initialize states
		self.states=False # Let this be a class of state variables and ode to evolve it.
		## Intialize other attributes
		for i in kwargs.iterkeys(): 
			vars(self)[i]=kwargs[i]

		## Initializing cell in space
		self.space=cellspace.space #pymunk space class
		inertia=pymunk.moment_for_poly(self.mass,\
		[(0,0), (0,self.height), (self.length,self.height), (self.length,0)], (0,0))
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
		self.space.add(self.body, self.box, self.left, self.right)
		##}}}

	def	_update_cell(self):
		##{{{
		"""
		Updates Inertia, Updates length, height etc.
		Does not update position.
		"""
		self.body.moment=pm.moment_for_poly(self.mass,[(0,0), (0,self.height), \
		(self.length,self.height), (self.height,0)], (0,0)) #Sets new inertia
		self.box.unsafe_set_vertices([(0,0), (0,self.height), (self.length,self.height), (self.length,0)],\
		(-self.length/2.,-self.height/2.))
		self.left.unsafe_set_offset((-self.length/2.,0))
		self.right.unsafe_set_offset((self.length/2.,0))
		self.left.unsafe_set_radius(self.radius)
		self.left.unsafe_set_radius(self.radius)
		##}}}
	
	def _update_states(self):
		##{{{
		"""
		Upadtes self.states. Stores them
		"""
		pass
		##}}}

	def _check_states(self):
		##{{{
		"""
		Depending on what the chemical state is, do sthg
		"""
		pass
		##}}}

	def _check_division(self):
		##{{{
		if self.length>self.division:
			self.divide()
		##}}}
	
	def evolve(self,dt):
		##{{{
		"""
		Helper function to decide what happens to the cell
		Move cell forward by one time step
		3 states. grow, quie, dead
		"""
		if self.cycle=='grow':
			self.grow(dt)
			self._check_division()
		elif self.cycle=='quie':
			pass
		elif self.cycle=='dead':
			pass
		##}}}

	def help(self):
		##{{{	
		"""
		To call for attributes from body, just use [self].body.position. Replace self with variable name
		"""
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
		self._update_cell()
		#Updating New Daughter
		self.cellspace.add_cell(pos2,self.body.angle,length=self.length) #add state next time

		##}}}

	def split_state(self):
		##{{{
		"""
		Splits the chemical state into 2. May add stochastics
		"""
		pass
		#Return a,b
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

	def get_pos(self):
		##{{{
		"""
		Gets the position of the cell
		"""
		pass
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

class cellspace():
	##{{{
	"""
	Keeps track of all cells. Decides which cells are active, which are not active. 
	"""
	def __init__(self,t0=0,dt=1.0/60.0,dim=(600,600),view=True):
		##{{{
		self.view=view # Determine if we launch pygame
		if self.view:
			pg.init()
			self.screen = pg.display.set_mode(dim)
			self.clock = pg.time.Clock()
		self.t0=t0 #Initial Time
		self.dt=dt #Time step
		self.time=self.t0
		self.space = pm.Space()
		self.cells={}
		self.active_cells={}
		self.inactive_cells={}
		self.count=0
		##}}}
		
	def add_cell(self,pos,angle=0,**kwargs):
		##{{{	
		"""
		Add cells to the space
		""" 
		num=self.count
		self.active_cells[num]=cell(self,num,pos,angle,**kwargs)
		self.count+=1
		##}}}

	def run(self,bg="black"):
		##{{{
		running=True
		counter=0
		while running:
			
			activecells=list(self.active_cells.itervalues())
			for cell in activecells:
				cell.evolve(self.dt)
				
			### Update physics
			self.space.step(self.dt)
			if np.mod(counter,100)==0:
				print "Number of cells=",len(self.active_cells)
	
			## Pygame
			if self.view:
				##Events
				for event in pg.event.get():
					if event.type == pg.QUIT:
						running = False
					elif event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE:
						running = False
					elif event.type == pg.KEYDOWN and event.key == pg.K_p:
						pg.image.save(screen, "contact_with_friction.png")
				## Draw on screen
				self.screen.fill(pg.color.THECOLORS[bg])
				pm.pygame_util.draw(self.screen, self.space)			
				self.clock.tick(50)
				pg.display.set_caption("fps: " + str(self.clock.get_fps()))
				pg.display.flip()
			## Misc
			self.time+=self.dt
			counter+=1
		##}}}
	##}}}

	

