import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_agg as agg

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

fig = plt.figure(figsize=[4, 4],dpi=100,frameon=False)
ax = fig.add_axes([0,0,1,1])
ax2 = fig.add_axes([0,0,1,1])
ax3 = fig.add_axes([0,0,1,1])
extent=(0,1,0,1)
img=ax.imshow(np.array([[1,0],[0,1]]),extent=extent,alpha=1,cmap=gen_cmap("#FFFF00"))
img2=ax2.imshow(np.array([[5,0,5],[0,5,0],[5,0,5]]),vmax=300,vmin=0,extent=extent,alpha=1*0.6,cmap=gen_cmap("#FF0000"))
img3=ax3.imshow(np.array([[1,0,1,0],[0,1,0,1],[1,0,1,0],[0,1,0,1]]),extent=extent,alpha=1*0.6*0.6,cmap=gen_cmap("#0000FF"))
ax.axis('off')
#ax.autoscale_view('tight')
#ax.set_frame_on(False)
#plt.tight_layout(pad=0)

canvas = agg.FigureCanvasAgg(fig)
canvas.draw()
renderer = canvas.get_renderer()
raw_data = renderer.tostring_rgb() #This is where information gets transferred to pygame
 
import pygame
from pygame.locals import *
 
pygame.init()
 
window = pygame.display.set_mode((400, 400), DOUBLEBUF)
screen = pygame.display.get_surface()
 
size = canvas.get_width_height()
 
surf = pygame.image.fromstring(raw_data, size, "RGB")
screen.blit(surf, (0,0))
pygame.display.flip()

#while True:
#	pass
