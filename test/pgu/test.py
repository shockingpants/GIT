import pygame
from pygame.locals import *
from pgu import gui
# The maximum frame-rate
FPS = 30
WIDTH,HEIGHT = 640,480


screen = pygame.display.set_mode((640,480),SWSURFACE)
class StarControl(gui.Table):
	def __init__(self,**params):
		gui.Table.__init__(self,**params)
		fg = (255,255,255)

		self.tr()
		self.td(gui.Label("Phil's Pygame GUI",color=fg),colspan=2)
		
		self.tr()
		self.td(gui.Label("Warp Speed: ",color=fg),align=1)
		self.td(gui.Switch(value=False,name='warp'))

def render(dt):
	warp = form['warp'].value #updates switch value
	if warp:
		print "It is selected."



#########################################
####		   Initialize
#########################################
form = gui.Form()
app = gui.App()
c = gui.Container(align=-1,valign=-1)
starCtrl = StarControl()
c.add(starCtrl,0,0)
app.init(c)
#######################################
clock = pygame.time.Clock()
done=False
while not done:
	for e in pygame.event.get():
		if e.type is QUIT: 
			done = True
		elif e.type is KEYDOWN and e.key == K_ESCAPE: 
			done = True
		else:
			app.event(e)

	# Clear the screen and render the stars
	dt = clock.tick(FPS)/1000.0
	screen.fill((0,0,0))
	render(dt)
	app.paint()
	pygame.display.flip()

