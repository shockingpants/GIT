#! /usr/bin/env python
import pygame
height=800
width=800
screen = pygame.display.set_mode((800, 800))
running=1
count=0
n=10
hint=height*1.0/n
wint=width*1.0/n
while running:
	event=pygame.event.poll()
	if event.type==pygame.QUIT:
		running = 0
	screen.fill((0,0,0))
	for i in xrange(n):
		pygame.draw.line(screen, (0, 0, 255), (0, 0+hint*i), (width-wint*i, 0))
	#pygame.draw.aaline(screen, (0, 255, 255), (639, 0), (0, 479))
	pygame.display.flip()
