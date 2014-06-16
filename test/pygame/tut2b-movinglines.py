#! /usr/bin/env python
# http://lorenzod8n.wordpress.com/category/pygame-tutorial/
import pygame
y = 0
dir = 1
running = 1
barheight = 124
screen = pygame.display.set_mode((800, 600));

barcolor = []
for i in range(1, 63):
    barcolor.append((0, 0, i*4))
for i in range(1, 63):
    barcolor.append((0, 0, 255 - i*4))

while running:
    event = pygame.event.poll()
    if event.type == pygame.QUIT:
        running = 0

    screen.fill((0, 0, 0))
    for i in range(0, barheight):
        pygame.draw.line(screen, barcolor[i], (0, y+i), (799, y+i))

    y += dir
    if y + barheight > 599 or y < 0:
        dir *= -1

    pygame.display.flip()
