###############################################################################################################
# CellSim
#
# Copyright Jonathan Teo   
# Contact   jonteo@mit.edu
# 2014
###############################################################################################################
import math, sys, random
import numpy as np
import pymunk as pm
import fipy as fp
import fipy.tools.numerix as fnumerix
#fipy, by national institute of standards and technology
from celltype import *
from viewer import *
import time

############################
####		Class
############################
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
		###### Viewer
		self.play=True #Variable specifying state
		self.pause=False
		if self.view:
			self.viewer=cellviewer(self) #This handles everything to do with viewing
			
		##}}}	
	def add_cell(self,celltyp,pos,angle=0,biochem=None,**kwargs):
		##{{{	
		"""
		Add cells to the space
		celltype should be the class itself (not an instance of the class)
		pos should be a tuple (x,y)
		Arguments specific to the cell type will be entered via kwargs
		TBD Need to change cell division in such a way that when the cell divides, certain new cell properties are maintained.
		""" 
		assert issubclass(celltyp, cell)
		num=self.count
		self.active_cells[num]=celltyp(self,num,pos,angle,biochem=biochem,**kwargs)
		self.count+=1
		return self.active_cells[num]
		##}}}
	def run(self,dt=None):
		"""
		if dt is none, use default dt
		"""
		##{{{
		counter=0
		dt=self.dt if dt is None else dt
		while self.play:
			if not self.pause:
				## Active Cell List	
				activecells=list(self.active_cells.itervalues())
				for cell in activecells:
					cell.evolve(self.dt)
					
				### Update physics
				self.space.step(self.dt)
				if np.mod(counter,100)==0:
					print "Number of cells=",len(self.active_cells)
		
				## Signaling Molecule
				self.solspace.run()
				
				## Misc
				self.time+=self.dt
				counter+=1
				
				## Pygame
				if self.view:
					self.viewer.plot()
			else:
				if self.view:
					self.viewer.plot()

			activecells=list(self.active_cells.itervalues())
			
			##}}}
	##}}}
class solspace(object):
##{{{
	"""
	A space where diffusable molecules diffuse using the various algorithm
	Advantages and Disadvantages --> http://www.cfm.brown.edu/people/sg/AM35odes.pdf
	Will use fipy as PDE solver
	"""
	def __init__(self,cellspac,dim=(100,100)):
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
		self.nx,self.ny=dim #Number of mesh along each axis
		self.dx=cs.dim[0]*cs.dx/1.0/self.nx #normalizes dx length to that of cellspace
		self.dy=cs.dim[1]*cs.dy/1.0/self.ny #Make sure it is float, not int
		#self.dx=dx=1.
		#self.dy=dy=1.
		Lx=self.nx*self.dx
		Ly=self.ny*self.dy
		self.mesh = fp.Grid2D(dx=self.dx, dy=self.dy, nx=self.nx, ny=self.ny)
		self.meshnum = self.mesh.getNumberOfCells() #Number of mesh grids in mesh
		self.species={} #dictionary of signal_species
		#self.viewer=_solviewer(self) #class responsible for all viewing
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
			oldpos[0].append(x*self.cellspace.dx)
			oldpos[1].append(y*self.cellspace.dy)
		#solpos=self.cellspace2solspace(oldpos)
		#ID=self.mesh._getNearestCellID(solpos)
		ID=self.mesh._getNearestCellID(oldpos)
		self.mask=np.zeros(self.meshnum)
		self.pos_ID=ID
		self.mask[ID]=1 #This is the actual mask of 1s and 0s, Can replace with True and False
		self.mask=self.mask>0
		return ID
		##}}}
	def _old_exchange(self):
		##{{{
		"""
		Goes through all the cells to get values and exchange it with solution
		solspace controls the exchange
		Assumption
		1) All cells exchange with the same kin and kout for a specific species of signaling mol
		Optimization
		1) same loop in _get_position and _exchange. Any way to combine them pythonically? 
		2) Too many loops, anyway to speed this up?
		"""
		dt=self.cellspace.dt
		#Get cell values
		int_val={} #Interior value
		kin={}
		kout={}
		for spec in self.species.itervalues():
			int_val[spec.name]=[]
		for cell in self.cellspace.active_cells.itervalues(): #Retrieves cells
			for spec in self.species.itervalues():
				int_val[spec.name].append(cell.biochem.get_latest(spec.name))	

		#Get solspace values
		ext_val={}
		assert "pos_ID" in vars(self) # Make sure to run _get_position before _exchange to generate ID
		for spec in self.species.itervalues():
			ext_val[spec.name]=spec.value[self.mask]

		#Actual exchange
		for spec in self.species.itervalues():
			kin=spec.kin
			kout=spec.kout
			dext_dt=kout*np.array(int_val[spec.name])-kin*np.array(ext_val[spec.name])
			int_val[spec.name]-=dext_dt*dt
			ext_val[spec.name]+=dext_dt*dt
			#Update
			spec.setValue(ext_val[spec.name], where=self.mask)
			for ind,cell in enumerate(self.cellspace.active_cells.itervalues()):
				cell.biochem.set_latest(spec.name,int_val[spec.name][ind])
		##}}}
	def _exchange(self):
		##{{{
		"""
		Goes through all the cells to get values and exchange it with solution
		solspace controls the exchange
		Assumption
		1) All cells exchange with the same kin and kout for a specific species of signaling mol
		Optimization
		1) same loop in _get_position and _exchange. Any way to combine them pythonically? 
		2) Too many loops, anyway to speed this up?
		"""
		dt=self.cellspace.dt

		for cell in self.cellspace.active_cells.itervalues(): #Retrieves cells
			oldpos=([],[])
			x,y=cell.get_pos()
			oldpos[0].append(x*self.cellspace.dx)
			oldpos[1].append(y*self.cellspace.dy)
			ID=self.mesh._getNearestCellID(oldpos)
			mask=np.zeros(self.meshnum)
			mask[ID]=1 #This is the actual mask of 1s and 0s, Can replace with True and False
			for spec in self.species.itervalues():
				#Get Value
				int_val=cell.biochem.get_latest(spec.name)
				ext_val=spec.value[ID]
				kin=spec.kin
				kout=spec.kout
				#Change
				dext_dt=kout*int_val-kin*ext_val
				int_val-=dext_dt*dt
				ext_val+=dext_dt*dt
				#Update VAlue
				cell.biochem.set_latest(spec.name,int_val)
				spec.setValue(ext_val, where=mask)
		##}}}
	def _update_interface(self):
		##{{{
		"""
		Changes the mesh CellVariable value based on position of cells and the rate of exchange.
		TDB: Write a more elegant solution that can be changed from a master script, not from the 
		__init__.py code.
		"""
		for spec in self.species.itervalues():
			spec.eq.solve(var=spec,dt=self.cellspace.dt)
			#spec.setValue(spec+0.2, where=self.mask) # This is just a test code.
		##}}}
	def add_species(self,name,degradation,diffusion,kin,kout,value=0.,**kwargs):
		##{{{
		"""
		Add class signal species
		# Note, it is a little different from cellspace.add_cell(). Here, the input 
		# is an instance of a class instead of referencing the class itself.

		-----------------------
		optional keywords = default --> description
		----------------------
		datamin = 0.0  --> lowest value represented on plot
		datamax = 1.0  --> highest value
		color = "0000FF" --> color representing that species.
		"""
		newspecies=signal_species(name,self.mesh,degradation,diffusion,kin,kout,value,**kwargs)
		newspecies.eq=fp.TransientTerm() == fp.DiffusionTerm(coeff=newspecies.diffusion)
		self.species[name]=newspecies #This a dictionary that keeps track of species.

		#datamin = kwargs['datamin'] if 'datamin' in kwargs else 0
		#datamax = kwargs['datamax'] if 'datamax' in kwargs else 1
		#if self.cellspace.view:
		#	self.viewer.add_species(spec=newspecies,datamin=datamin, datamax=datamax) #TEMP
		
		##}}}
	def run(self):
		##{{{
		"""
		Evolve the system
		"""
		#self._get_position()
		self._update_interface()
		self._exchange()

		
		##}}}
	def cellspace2solspace(self,pos):
		##{{{
		"""
		NOTE: USELESS RIGHT NOW
		converts cellspace coordinate into solspace coordinate
		pos is a tuple of (x,y)
		"""
		oldx=np.array(pos[0])*self.cellspace.dx
		oldy=np.array(pos[1])*self.cellspace.dy
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
	value is the initial value
	"""
	def __init__(self, name, mesh, degradation ,diffusion, kin, kout, value=0., **kwargs):
		##{{{
		super(signal_species, self).__init__(name=name, mesh=mesh, value=value)
		self.degradation=degradation
		self.diffusion=diffusion
		self.name=name
		self.kwargs=kwargs
		self.kin=kin
		self.kout=kout
		##}}}
	##}}}
class GUI(object):
	##{{{
	"""
	pass
	"""
	def __init__():
		pass
	##}}}
#### Need to configure space for constraints, but work with no constraints for now. 
##}}}

#############################
####  Utility Function
#############################


