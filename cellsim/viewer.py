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

import matplotlib
matplotlib.use("Agg")
import matplotlib.backends.backend_agg as agg
import matplotlib.pyplot as plt

import fipy as fp
import fipy.tools.numerix as fnumerix

from pgu import gui

class cellviewer(object):
	##{{{
	"""
	Creates a plot object that plots directly onto pygame
	"""
	def __init__(self,cs):
		##{{{
		"""
		>>>type(cs) 
		cellspace
		"""
		self.cellspace=cs
		self.solspace=sols=cs.solspace
		# Initialize pygame
		pg.init()
		self.screen = pg.display.set_mode(cs.dim)
		self.clock = pg.time.Clock()
		# Initialize MPL onto pygame
		dpi=100.
		self.fig=plt.figure(figsize=[cs.dim[0]/dpi,cs.dim[1]/dpi],dpi=dpi,frameon=False)
		self.ax=self.fig.add_axes([0,0,1,1])
		self.ax.axis("off")
		self.extent = (0,sols.nx*sols.dx,0,sols.ny*sols.dy) #xmin,ymin,xmax,ymax
		self.canvas = agg.FigureCanvasAgg(self.fig)
		self.size = self.canvas.get_width_height() # Needed as an argument for pygame.image.fromstring
		self.gui = Control(cs,cs.solspace)
		self.init_state=False
		#self.img=self.ax.imshow(np.zeros(2),extent=self.extent,vmin=datamin,vmax=datamax) #Just for initializing. Zero data
		##}}}
	def _initialize(self):
		##{{{
		"""
		Initializes after signaling species have been added
		"""
		##------ Signaling Species Specific
		self.img={} #Dictionary referencing the plot
		for ind,spec in enumerate(self.solspace.species.itervalues()):
			data=fnumerix.reshape(fnumerix.array(spec), spec.mesh.shape[::-1])[::-1]
			kwargs=spec.kwargs
			datamin = kwargs['datamin'] if 'datamin' in kwargs else 0
			datamax = kwargs['datamax'] if 'datamax' in kwargs else 1
			color = kwargs['color'] if 'color' in kwargs else "#0000FF"
			self.img[spec.name]=self.ax.imshow(data,extent=self.extent,alpha=1*0.5**ind,vmin=datamin,vmax=datamax,cmap=gen_cmap(color)) 
			self.gui.add_switch(spec.name) # Add switch to GUI
		self.gui._initialize()
		self.init_state=True
		##}}}
	def _plot_cells(self):
		##{{{
		"""
		Plots the layer of cells on pygame
		"""
		for event in pg.event.get():
			if event.type == pg.QUIT:
				self.cellspace.play = False
			elif event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE:
				self.cellspace.play = False
			elif event.type == pg.KEYDOWN and event.key == pg.K_p:
				pg.image.save(self.screen, "contact_with_friction.png")
			else:
				self.gui.app.event(event) #Sends event into gui
		## Draw on screen
		self.screen.fill(pg.color.THECOLORS["white"])
		activecells=list(self.cellspace.active_cells.itervalues())
		for cell in activecells:
			cell.draw()
		self.clock.tick(50)
		pg.display.set_caption("fps: " + str(self.clock.get_fps()))
		##}}}
	def _plot_sol(self):
		##{{{
		"""
		Plots the solution layer
		plots based on GUI ticks		
		"""
		#Update
		for spec in self.solspace.species.itervalues():
			if self.gui.form[spec.name].value: #Checks the checkbox/gui/switch status
				data=fnumerix.reshape(fnumerix.array(spec), spec.mesh.shape[::-1])[::-1]
				self.img[spec.name].set_data(data)
			else:
				self.img[spec.name].set_data(np.array([[0],[0]]))
			
		self.canvas.draw()
		renderer = self.canvas.get_renderer()
		raw_data = renderer.tostring_rgb()
		surf = pg.image.fromstring(raw_data, self.size, "RGB")
		self.screen.blit(surf, (0,0))
		#self.fig.canvas.draw()
		##}}}
	def plot(self):
		##{{{
		"""
		TBD: Manage ticks and GUI
		"""
		if not self.init_state:
			self._initialize()
		self._plot_cells()
<<<<<<< HEAD
		if len(self.solspace.species)>0:
			self._plot_sol()		
=======
		boolean=[self.gui.form[spec.name].value for spec in self.solspace.species.itervalues()] #list of switch status
		if len(self.solspace.species)>0 and any(boolean):
			self._plot_sol()
		self.gui.app.paint()
>>>>>>> FETCH_HEAD
		pg.display.flip()
		##}}}

	##}}}
class Control(gui.Table):
	##{{{
	"""
	GUI coded using pgu.
	"""
	def __init__(self,cellspace,solspace,**params):
		##{{{
		gui.Table.__init__(self,**params)
		self.cellspace=cellspace
		self.solspace=solspace
		self.fg = (0,0,0)
		self.add_play()
		self.add_pause()
		self.tr()
		self.td(gui.Label("Signal layers",color=self.fg),colspan=2)
		self.form = gui.Form()
		self.app = gui.App()
		self.c = gui.Container(align=-1,valign=-1)

		##}}}
	def _initialize(self):
		##{{{
		"""
		This is different from init. Adds switches once signal species have been specified
		"""
		self.c.add(self,0,0)
		self.app.init(self.c)

		##}}}
	def add_switch(self,label,value=False):
		##{{{
		"""
		Adds a check box on which signal layer to display
		"""
		self.tr()
		self.td(gui.Label("{0:s}: ".format(label),color=self.fg),align=1)
		self.td(gui.Switch(value=value,name=label))
		##}}}
	def add_play(self):
		##{{{
		"""
		Adds a play button
		"""
		but = gui.Button("Play")
		but.connect(gui.CLICK,self.set_pause,False)
		self.tr()
		self.td(but)
		##}}}

	def add_pause(self):
		##{{{
		"""
		Adds a pause button
		"""
		but = gui.Button("Pause")
		but.connect(gui.CLICK,self.set_pause,True)
		self.tr()
		self.td(but)
		##}}}
		
	def set_pause(self,state):
		"""
		state --> True or False
		Edits the pause variable in cellspace
		"""
		self.cellspace.pause=state
	##}}}
##############################
####  Utility Function
#############################
def hex_to_rgb(value):
	##{{{
	"""
	http://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa
	"""
	value = value.lstrip('#')
	lv = len(value)
	return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))
	##}}}
def gen_cmap(color):
	##{{{
	"""
	Takes in tuple of rgb value or hex value
	Example:	
	>>> gen_cmap(#FF0000) 
	>>> gen_cmap((255,0,0))		
	Returns:
	cdict	 #A matplotlib color construct
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
	RGB=[i/256. for i in RGB]#Normalize rgb
	cdict = {'red':   ((0.0,  1.0, 1.0),\
			(1.0,  RGB[0], RGB[0])),\
			'green': ((0.0,  1.0, 1.0),\
			(1.0,  RGB[1], RGB[1])),\
			'blue':  ((0.0,  1.0, 1.0),\
			(1.0,  RGB[2], RGB[2]))}
	return matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	##}}}

