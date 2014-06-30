import math, sys, random
import pygame as pg
import pygame.locals as pl
import pygame.color as pc
import pymunk as pm
from pymunk import Vec2d
import pygame_draw
import numpy as np
import fipy as fp
############################
####		Class
############################
##}}}
class cell(object):
	##{{{
	"""
	Every attribute of a cell will be written here.
	Should not be called directly since it contains no shape information. Should call using subclasses
	"""
	def __init__(self,cellspace,num,pos,angle,**kwargs):
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
		self.states=False # Let this be a class of state variables and ode to evolve it.
		## Intialize other attributes
		for i in kwargs.iterkeys(): 
			vars(self)[i]=kwargs[i]
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
class cellspace(object):
	##{{{
	"""
	Keeps track of all cells. Decides which cells are active, which are not active. 
	2nd space is the space that pymunk use to keep track of the body and shapes (physics space).
	3rd space is the space that screen that pygame use to draw shapes.
	"""
	def __init__(self,t0=0,dt=1.0/60.0,dim=(600,600),view=True):
		##{{{
		"""
		each pixel is 0.1um
		"""
		self.view=view # Determine if we launch pygame
		self.dx=self.dy=0.1 #0.1um
		if self.view:
			pg.init()
			self.screen = pg.display.set_mode(dim)
			self.clock = pg.time.Clock()
		self.dim=dim
		self.t0=t0 #Initial Time
		self.dt=dt #Time step
		self.time=self.t0
		self.space = pm.Space()
		self.cells={}
		self.active_cells={}
		self.inactive_cells={}
		self.count=0
		self.solspace=solspace(self)
		##}}}	
	def add_cell(self,celltype,pos,angle=0,**kwargs):
		##{{{	
		"""
		Add cells to the space
		celltype should be the class itself (not an instance of the class)
		pos should be a tuple (x,y)
		Arguments specific to the cell type will be entered via kwargs
		TBD Need to change cell division in such a way that when the cell divides, certain new cell properties are maintained.
		""" 
		assert issubclass(celltype, cell)
		num=self.count
		self.active_cells[num]=celltype(self,num,pos,angle,**kwargs)
		self.count+=1
		return self.active_cells[num]
		##}}}
	def run(self,bg="white"):
		##{{{
		running=True
		counter=0
		while running:
			## Active Cell List	
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
						pg.image.save(self.screen, "contact_with_friction.png")
				## Draw on screen
				self.screen.fill(pg.color.THECOLORS[bg])
				for cell in activecells:
					cell.draw()
				#pygame_draw.draw(self.screen, self.space)			
				self.clock.tick(50)
				pg.display.set_caption("fps: " + str(self.clock.get_fps()))
				pg.display.flip()

			## Signaling Molecule
			self.solspace.run()

			## Misc
			self.time+=self.dt
			counter+=1
		##}}}
	##}}}
class solspace(object):
##{{{
	"""
	A space where diffusable molecules diffuse using the various algorithm
	Advantages and Disadvantages --> http://www.cfm.brown.edu/people/sg/AM35odes.pdf
	Will use fipy as PDE solver
	"""
	def __init__(self,cellspac,dim=(200,200)):
		##{{{
		"""
		cellspac is an instance of cellspace that this solspace will be associated with
		We dont want the PDE solver to be huge, so the mesh size should be small, probably about 200x200
		There may be other considerations such as grid size to cell size ratio, but for now, assume 200x200
		Dimension of solspace will typically not match dimension of cellspace
		"""
		assert isinstance(cellspac,cellspace)
		#------Mesh
		cs=self.cellspace=cellspac #refers to the cell space that contains this sol space
		self.dim=dim
		nx,ny=dim #Number of mesh along each axis
		#dx=nx/1.0/cs.dim[0]*cs.dx #normalizes dx length to that of cellspace
		#dy=ny/1.0/cs.dim[1]*cs.dy
		self.dx=dx=1.
		self.dy=dy=1.
		Lx=nx*dx
		Ly=ny*dy
		self.mesh = fp.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
		self.meshnum = self.mesh.getNumberOfCells() #Number of mesh grids in mesh
		self.species={} #dictionary of signal_species
		self.viewer={} #TEMP dictionary of viewers
		##}}}
	def _update_interface(self):
		##{{{
		"""
		Changes the mesh CellVariable value based on position of cells and the rate of exchange.
		TDB: Write a more elegant solution that can be changed from a master script, not from the 
		__init__.py code.
		"""
		for spec in self.species.itervalues():
			print spec.eq
			print self.dim
			print self.dx
			print self.cellspace.dt
			raw_input()
			spec.eq.solve(var=spec,dt=self.cellspace.dt)
			spec.setValue(spec+1, where=self.mask) # This is just a test code.
			self.viewer[spec.name].plot()
		##}}}
	def _get_position(self):
		##{{{
		"""
		Uses position of cells.
		Returns mask, input for fp.CellVariable.setValue
		"""
		oldpos=([],[])
		for cell in self.cellspace.active_cells.itervalues(): #Retrieves cells
			x,y=cell.get_pos()
			oldpos[0].append(x)
			oldpos[1].append(y)
		solpos=self.cellspace2solspace(oldpos)
		ID=self.mesh._getNearestCellID(solpos)
		self.mask=np.zeros(self.meshnum)
		self.mask[ID]=1 #This is the actual mask of 1s and 0s, Can replace with True and False	
		##}}}
	def add_species(self,name,degradation,diffusion,value=0.,**kwargs):
		##{{{
		"""
		Add class signal species
		# Note, it is a little different from cellspace.add_cell(). Here, the input 
		# is an instance of a class instead of referencing the class itself.
		"""
		newspecies=signal_species(name,self.mesh,degradation,diffusion,value,**kwargs)
		newspecies.eq=fp.TransientTerm() == fp.DiffusionTerm(coeff=newspecies.diffusion)
		self.species[name]=newspecies #This a dictionary that keeps track of species.
		self.viewer[name]=fp.Viewer(vars=newspecies, datamin=0., datamax=1) #TEMP
		
		##}}}
	def run(self):
		##{{{
		"""
		Evolve the system
		"""
		self._get_position()
		self._update_interface()

		
		##}}}
	def cellspace2solspace(self,pos):
		##{{{
		"""
		converts cellspace coordinate into solspace coordinate
		pos is a tuple of (x,y)
		"""
		oldx=np.array(pos[0])
		oldy=np.array(pos[1])
		cs=self.cellspace
		x=oldx/cs.dim[0]*self.dim[0]
		y=oldy/cs.dim[1]*self.dim[1]
		return (x,y)
		##}}}
	def solspace2cellspace(self,pos):
		##{{{
		"""
		converts solspace coordinate into cellspace coordinate
		pos is a tuple of (x,y)
		"""
		cs=self.cellspace
		x=pos[0]/self.dim[0]*cs.dim[0]
		y=pos[1]/self.dim[1]*cs.dim[1]
		return (x,y)
		##}}}
	
	
##}}}
class signal_species(fp.CellVariable):
	##{{{
	"""
	Create a new subclass using fp.CellVariable
	"""
	def __init__(self, name, mesh, degradation ,diffusion, value=0., **kwargs):
		super(signal_species, self).__init__(name=name, mesh=mesh, value=value)
		self.degradation=degradation
		self.diffusion=diffusion
		self.name=name
		self.kwargs=kwargs
	##}}}
##--------Custom Cells--------------	
class cellp(cell):
	##{{{
	"""
	Cell approximated by a polygon
	1 px is 0.1um	
	"""
	def __init__(self,cellspace,num,pos,angle,radius=4,length=8,vertices=6,mass=0.1,**kwargs):
		##{{{
		## Cell Attributes
		super(cellp, self).__init__(cellspace,num,pos,angle,**kwargs)
		self.radius=radius #Radius of hemisphere at cell end
		self.length=length #The length of the cell measured from the center of the two hemisphere
		self.height=self.radius*2
		self.vertices=vertices #Number of vertices in the hemisphere
		self.mass=mass #Self Explanatory
		self.color=pg.color.THECOLORS["red"]
		self.growthrate=0.05 #Growth rate can be a function of a state variable eventually
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
		Updates Inertia, Updates length, height etc.
		Does not update position.
		"""
		self._genpoly()
		self.body.moment=pm.moment_for_poly(self.mass,self.ver, (0,0)) #Sets new inertia
		self.box.unsafe_set_vertices(self.ver,(-self.length/2.,-self.height/2.))
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
		# Change to grow, quienescence, death
		# Change from active list to inactive list
		pass
		##}}}
	def _check_division(self):
		##{{{
		self.division=5*self.radius #Threshold for division
		if self.length>self.division:
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
		screen=self.cellspace.screen
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

		###### Gets position of daughter cells
		oldpos=self.body.position
		a=self.box.get_vertices()[5] #top left hand corner
		b=self.box.get_vertices()[13]
		c=self.box.get_vertices()[12]
		d=self.box.get_vertices()[6]
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
		temp=self.cellspace.add_cell(cellp,pos2,self.body.angle,length=self.length,radius=self.radius,**self.kwargs) #TBDadd state next time
		temp.body.velocity=self.body.velocity
		temp.body.angular_velocity=self.body.angular_velocity
		# TAKENOTE, MAKE SURE THERE ARE NO NEW KEYWORDS ADDED, OTHERWISE **KWARGS WILL INCREASE WITH EVERY DIVISION
		# AND WE DON'T WANT THAT

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
		# Change to grow, quienescence, death
		# Change from active list to inactive list
		pass
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
		screen=self.cellspace.screen
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
		"""
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
		self._update_cell()
		#Updating New Daughter
		self.cellspace.add_cell(cellc,pos2,self.body.angle,length=self.length,radius=self.radius) #add state next time

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

	##}}}
##-------Custop cellspace-------------	
#### Need to configure space for constraints, but work with no constraints for now. 
##}}}

#############################
####  Utility Function
#############################
##{{{
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
##}}}
	

