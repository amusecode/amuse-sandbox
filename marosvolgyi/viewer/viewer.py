try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    import pygame
    from pygame.locals import *
except ImportError:
    if __name__ == '__main__':
        raise Exception("OpenGL and pygame bindings are not installed, cannot run")
    
        
#from OpenGL.GLUT import *


from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system 
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.kingmodel import MakeKingModel
from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.octgrav.interface import Octgrav

from datetime import date, timedelta
import time

import sys
import numpy as np
import random

class Viewer(object):
    def __init__(self, *kargs):
        if len(kargs)==2:
            display_size = kargs[1]
            self.stars = kargs[0]
        elif len(kargs)==1:
            display_size = kargs[0]
            self.stars = []

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
            #gluSphere(pObj, r[0], 52, 26)# Bottom
            gluSphere(pObj, r[0], 8, 6)# Bottom

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
                    pygame.image.save(self.screen, 'screenshot.png')
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
        R = [[0.0005,0,0,0]]

        for i, s in enumerate(stars):
            r = s.position
            x = r[0].value_in(nbody_system.length)
            y = r[1].value_in(nbody_system.length)
            z = r[2].value_in(nbody_system.length)
            R.append([0.01, x, y, z])

        self.draw_scene(R)
        self.draw_plane()
        pygame.display.flip()

    def animate(self):

        if self.autoevolve:
            self.days += 5
            
        R = [[0.1, 0,0,0]]

        for i in self.stars:
            r, v = i.get_vectors_at_date(self.days)
            R.append([0.05, r[0], r[1], r[2]])

        self.omega += 0.05
        self.render(R)

if __name__ == "__main__":

    nstars = int(sys.argv[1])
    workers = int(sys.argv[2])
    method = sys.argv[3]
    print nstars
    seed = None
    #convert_nbody = nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.AU)
    stars = MakePlummerModel(nstars).result#, convert_nbody, random_state = seed).result
    #stars2= MakePlummerModel(nstars).result
    #stars2 = MakeKingModel(nstars, None, W0=6).result
    #stars2.x += 0.5 | nbody_system.length
    #stars2.vy = 35 | nbody_system.speed
    #stars2.mass = 1.0 |nbody_system.mass
    #stars.add_particles(stars2)
    s = Viewer(stars ,(1024, 780))

    #stars.scale_to_standard()
    if method == 'octgrav':
        gravity = Octgrav()
    elif method == 'phigrape':					
        gravity = PhiGRAPE(PhiGRAPE.NBODY)
    elif method == 'bhtree':
        gravity = BHTree(number_of_workes = workers) 
    elif method == 'hermite':
        gravity = Hermite(
                          number_of_workers = workers
                          #debugger = "xterm",
                          #redirection = "none"
                          )
    gravity.initialize_code()
    gravity.parameters.epsilon_squared = 0.001 | nbody_system.length ** 2
    
    stars.radius = 0.000001 | nbody_system.length
    gravity.particles.add_particles(stars)
    #gravity.stopping_conditions.pair_detection.enable()
    from_model_to_gravity = stars.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(stars)
    
    gravity.commit_particles()
    model_time = 0.0
    pairs = 0
    while s.go:
        flag = time.time()
        #s.animate()
        model_time +=0.01

        gravity.evolve_model(model_time|nbody_system.time)
        from_gravity_to_model.copy()
        rterr = 0
        while (time.time() - flag)<0.01:
            rterr +=1
        if rterr == 1: 
            print "realtime error"

        #if gravity.stopping_conditions.pair_detection.is_set():
        #    pairs +=1
        #    print "Pair % i" % pairs
        #    print len(gravity.stopping_conditions.pair_detection.particles(0).key)

        s.renderamuse(stars)
        s.handle_events()
        
    gravity.stop()
