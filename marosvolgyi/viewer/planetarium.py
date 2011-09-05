try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    import pygame
    from pygame.locals import *
except ImportError:
    if __name__ == '__main__':
        raise Exception("OpenGL and pygame bindings are not installed, cannot run")
    
from amuse.io.horizons import LoadStar
from amuse.units import units
from datetime import date, timedelta
import time

import sys
import numpy as np

from amuse.support import data
class SolarSystemView(object):
    def __init__(self, *kargs):
        if len(kargs)==2:
            list_of_planets = kargs[0] 
            display_size = kargs[1]
            self.planets = list_of_planets
        elif len(kargs)==1:
            display_size = kargs[0]
            self.planets = []

        pygame.init()

        self.screen = pygame.display.set_mode(display_size,  HWSURFACE|OPENGL|DOUBLEBUF)
        #if you leave the last argument you get a segmentation fault!
        pygame.display.set_caption('Solar system viewer')

        self.gl_init(display_size)

        pygame.key.set_repeat(10, 5)#keyboard repeat delay, repeat in ms

        self.drag = False
        self.trail = False
        self.box = False
        self.autoevolve = False

        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.omega = 0.0
        self.zoom = -2
        self.days = 0
        self.Trans = [0,0]
        self.go = True

    def add_planet(self, planet):
        self.planets.append(planet)

    def remove_panet(self, planet_name):
        index = [i for i, planet in enumerate(self.planets) if planet.name==planet_name]
        if len(index)==0:
            print "No such planet, did not remove anything"
        else:
            self.planets.pop(index)

    def gl_init(self, display_size):

        glViewport(0, 0, display_size[0], display_size[1])
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60.0, 1.0*display_size[0]/display_size[1], .1, 1000.)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        glEnable(GL_DEPTH_TEST)

        glShadeModel(GL_FLAT)
        glClearColor(.0, .0, .0, 0.0)

        glEnable(GL_COLOR_MATERIAL)

        glEnable(GL_DEPTH_TEST)# Hidden surface removal
        glFrontFace(GL_CCW)# Counterclockwise polygons face out
        glEnable(GL_CULL_FACE) 

        glEnable(GL_LIGHTING)
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  (1, 1, 1, 1))            
        glLightfv(GL_LIGHT0, GL_POSITION,  (0, 1, 15, 1))
        #glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        #glMaterialfv(GL_FRONT, GL_SPECULAR,(1,1,1,1))
        #glMateriali(GL_FRONT,GL_SHININESS,128)

    def draw_plane(self):
        glPushMatrix()
        glEnable(GL_DEPTH_TEST)
        pObj = gluNewQuadric()
        gluQuadricNormals(pObj, GLU_LINE)
        gluQuadricDrawStyle(pObj, GLU_LINE)
        #Main Body
        glColor3f(.0, .0, 1)
        glColor
        glLoadIdentity()
        glTranslate(self.Trans[0],self.Trans[1],self.zoom)
        glRotatef(self.alpha, 1.0, 0.0, 0.0)
        glRotatef(self.beta, 0.0, 1.0, 0.0)
        glRotatef(self.gamma, 0.0, 0.0, 1.0)
        gluDisk(pObj, 0.0, 50.0, 50,30)
        glPopMatrix()

    def draw_scene(self, R):
        #Clear the window with current clearing color
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        #Move object back and do in place rotation

        #Draw something
        glPushMatrix()

        for i, r in enumerate(R):
            glEnable(GL_DEPTH_TEST)
            pObj = gluNewQuadric()
            gluQuadricNormals(pObj, GLU_LINE)
            gluQuadricDrawStyle(pObj, GLU_LINE)
            #Main Body
            glColor3f(1.0, 1.0, 1.0)
            glLoadIdentity()

            glTranslate(self.Trans[0],self.Trans[1],self.zoom)
            glRotatef(self.alpha, 1.0, 0.0, 0.0)
            glRotatef(self.beta, 0.0, 1.0, 0.0)
            glRotatef(self.gamma, 0.0, 0.0, 1.0)
            glTranslate(r[1], r[2],r[3])
            gluSphere(pObj, r[0], 52, 26)# Bottom
            
        glPopMatrix()
        glFlush()

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
            elif event.type == pygame.MOUSEMOTION and self.drag:
                dx, dy = event.rel
                if dy<0:
                    self.alpha-=1
                if dy>0:
                    self.alpha+=1
                if dx<0:
                    self.beta-=1
                if dx>0:
                    self.beta+=1
                
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 4:
                    self.zoom +=0.5
                elif event.button == 5:
                    self.zoom -=0.5
                else:
                    self.drag = True
            elif event.type == pygame.MOUSEBUTTONUP:
                self.drag = False
            elif event.type == KEYDOWN:
                if event.key == K_q:
                    self.go = False
                if event.key == K_UP:
                    self.Trans[1] -= .1
                if event.key == K_DOWN:
                    self.Trans[1] += .1
                if event.key == K_LEFT:
                    self.Trans[0] += .1
                if event.key == K_RIGHT:
                    self.Trans[0] -= .1
                if event.key == K_KP8:
                    self.alpha += 1
                if event.key == K_KP4:
                    self.beta -= 1
                if event.key == K_KP6:
                    self.beta += 1
                if event.key == K_KP2:
                    self.alpha -= 1
                if event.key == K_KP7:
                    self.gamma -= 1
                if event.key == K_KP9:
                    self.gamma += 1
                if event.key == K_KP_PLUS:
                    self.zoom += 0.3
                if event.key == K_KP_MINUS:
                    self.zoom -= 0.3
                if event.key == K_a:
                    self.days -= 10
                    if self.days <0:
                        self.days = 0
                if event.key == K_l:
                    self.days += 1
                    if self.days > 60000:
                        self.days = 0
                if event.key == K_s:
                    #screenshot
                    pygame.image.save(screen, 'screenshot.png')
                if event.key == K_e:
                    print "Auto evolve mode toggle"
                    self.autoevolve = not self.autoevolve
                    print self.autoevolve
                if event.key == K_t:
                    print "Toggle trails"
                    self.trail = not self.trail
                    print self.trail
                if event.key == K_d:
                    print date(1971,10,26)+timedelta(self.days)
                if event.key == K_h:
                    print "Help info...none yet :("

    def render(self, R):
        self.draw_scene(R)
        self.draw_plane()
        pygame.display.flip()

    def renderamuse(self, stars):
        R = [[0.05,0,0,0]]

        for i, s in enumerate(stars):
            r = s.position
            x = r[0].value_in(units.AU)
            y = r[1].value_in(units.AU)
            z = r[2].value_in(units.AU)
            R.append([0.1, x, y, z])

        self.draw_scene(R)
        self.draw_plane()
        pygame.display.flip()

    def animate(self):

        if self.autoevolve:
            self.days += 5
            
        R = [[0.1, 0,0,0]]

        for i in self.planets:
            r, v = i.get_vectors_at_date(date(1974,10,26)+timedelta(self.days))
            R.append([0.05, r[0], r[1], r[2]])

        self.omega += 0.05
        self.render(R)

if __name__ == "__main__":
    
    list_of_planets = []

    for i in ['earth', 'moon','mars','jupiter','saturn','neptune','uranus','pioneer10','pioneer11']:
        try:
            list_of_planets.append(LoadStar(i))
        except:
            print "No such planet found, check if horizons files exist or if the names are spelled correctly"
            sys.exit()

    S = SolarSystemView(list_of_planets,(1200,800))

    while S.go:
        flag = time.time()
        while time.time()-flag<0.01:
            S.animate()
        S.handle_events()
        

