import sys

import pygame
from pygame.locals import *
from pygame.color import *

import pymunk
from pymunk.vec2d import Vec2d
from pymunk.pygame_util import draw, from_pygame

width = height = 600
def main():
    ### PyGame init
    pygame.init()
    screen = pygame.display.set_mode((width,height))
    clock = pygame.time.Clock()
    running = True
    font = pygame.font.SysFont("Arial", 16)

    ### Physics stuff
    space = pymunk.Space()

    # walls - the left-top-right-bottom walls
    static_lines = [pymunk.Segment(space.static_body, (50, 50), (50, 550), 5)
                ,pymunk.Segment(space.static_body, (50, 550), (550, 550), 5)
                ,pymunk.Segment(space.static_body, (550, 550), (550, 50), 5)
                ,pymunk.Segment(space.static_body, (50, 50), (550, 50), 5)
                ]
    space.add_static(static_lines)

    ### "Cannon" that can fire arrows
    cannon_body = pymunk.Body()
    player_shape = pymunk.Circle(cannon_body, 25)
    cannon_body.position = 100,100

    space.add(player_shape)

    while running:
        for event in pygame.event.get():
            if event.type == QUIT or \
                event.type == KEYDOWN and (event.key in [K_ESCAPE, K_q]):
                running = False

        mpos = pygame.mouse.get_pos()
        p = from_pygame( Vec2d(mpos), screen )
        mouse_position = p
        cannon_body.angle = (mouse_position - cannon_body.position).angle

        ### Clear screen
        screen.fill(pygame.color.THECOLORS["black"])

        ### Draw stuff
        draw(screen, space)

        ### Update physics
        fps = 60
        dt = 1./fps
        space.step(dt)

        ### Info and flip screen
        screen.blit(font.render("fps: " + str(clock.get_fps()), 1, THECOLORS["white"]), (0,0))
        screen.blit(font.render("Aim with mouse, click to fire", 1, THECOLORS["darkgrey"]), (5,height - 35))
        screen.blit(font.render("Press R to reset, ESC or Q to quit", 1, THECOLORS["darkgrey"]), (5,height - 20))

        pygame.display.flip()
        clock.tick(fps)

if __name__ == '__main__':
    sys.exit(main())
