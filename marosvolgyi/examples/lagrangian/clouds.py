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
from amuse.community.hermite0.interface import Hermite
from amuse.community.ph4.interface import ph4
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.octgrav.interface import Octgrav

from datetime import date, timedelta
import time

import sys
import numpy as np
import random

class Viewer(object):
    def __init__(self, *kargs):
        if len(kargs)==3:
            display_size = kargs[1]
            self.stars = kargs[0]
            self.converter = kargs[2]
        if len(kargs)==2:
            display_size = kargs[1]
            self.stars = kargs[0]
        elif len(kargs)==1:
            display_size = kargs[0]
            self.stars = []

        self.speed = 0
        self.dx = 0
        self.dy=0
        self.satposx = 0
        self.satposy = 0
        self.cos = False

        pygame.init()

        self.screen = pygame.display.set_mode(display_size,  HWSURFACE|OPENGL|DOUBLEBUF)
        #if you leave the last argument you get a segmentation fault!
        pygame.display.set_caption('Solar system viewer')

        self.gl_init(display_size)

        pygame.key.set_repeat(1000, 5)#keyboard repeat delay, repeat in ms

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

        glTranslate(self.Trans[0], self.Trans[1] ,self.zoom)

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

            glTranslate(self.Trans[0], self.Trans[1],self.zoom)
            glRotatef(self.alpha, 1.0, 0.0, 0.0)
            glRotatef(self.beta, 0.0, 1.0, 0.0)
            glRotatef(self.gamma, 0.0, 0.0, 1.0)
            glTranslate(r[1], r[2],r[3])
            #gluSphere(pObj, r[0], 52, 26)# Bottom
            gluSphere(pObj, r[0], 12, 10)# Bottom

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
                    self.zoom += 0.02
                if event.key == K_KP_MINUS:
                    self.zoom -= 0.02

                if event.key == K_p:
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
                if event.key == K_h:
                    print "Help info...none yet :("

    def render(self, R):
        self.draw_scene(R)
        self.draw_plane()
        pygame.display.flip()

    def renderamuse(self, stars, filename):
        R = [[.1, 0, 0, 0]]

        for i, s in enumerate(stars):
            r = s.position
            
            x = r[0].value_in(units.AU)
            y = r[1].value_in(units.AU)
            z = r[2].value_in(units.AU)
            R.append([s.radius.value_in(units.km), x, y, z])

        self.draw_scene(R)
        self.draw_plane()
        pygame.display.flip()
        #pygame.image.save(self.screen, filename)

    def animate(self):

        if self.autoevolve:
            self.days += 5
            
        R = [[0.1, 0,0,0]]

        for i in self.stars:
            r, v = i.get_vectors_at_date(self.days)
            R.append([0.05, r[0], r[1], r[2]])

        self.omega += 0.05
        self.render(R)

def diffangle(v,w):
    vn = np.dot(v,v.transpose())[0,0]**0.5
    wn = np.dot(w,w.transpose())[0,0]**0.5
    a = np.arctan2(v[0,1],v[0,0])-np.arctan2(w[0,1],w[0,0])
    return a

def rotate(vector, phi, teta, psi):
    
    normed_vect = [(i).value_in(units.AU) for i in vector]
    vt = np.matrix(normed_vect).transpose()
    
    Rotx = np.matrix([[1, 0, 0],[0,np.cos(phi), -np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
    Roty = np.matrix([[np.cos(teta), 0, np.sin(teta)],[0,1, 0],[-np.sin(teta),0,np.cos(teta)]])
    Rotz = np.matrix([[np.cos(psi), -np.sin(psi), 0],[np.sin(psi),np.cos(psi), 0],[0,0,1]])
    Rot = Rotx*Roty*Rotz

    vtt = 1.447515189957170E-02 * np.array((Rot * vt).transpose())
    return units.AUd(vtt)

def remove_outershell(dust):
    outer = dust.select(lambda r: r.norm()>2.0 | units.AU,["position"])
    dust.remove_particles(outer)

def generate_dust(dust):
    #assuming two AU sphere
    n_particles = len(dust)
    initial_pos = 5.0*(np.random.random([n_particles, 3])-0.5) | units.AU
    initial_pos[::,2] *=0.01
    dust.position = initial_pos
    dust.velocity = rotate(dust.position, 0, 0, 0.5*3.1415)
    
    #dust.velocity = np.zeros([n_particles, 3])|units.AUd
    dust.mass = 1000*np.ones(n_particles)|units.kg
    dust.radius = 0.01*np.ones(n_particles)|units.km
    remove_outershell(dust)
    
def remove_core_dust(dust):
    inner = dust.select(lambda r: r.norm()<0.5 | units.AU,["position"])
    dust.remove_particles(inner)

def remove_fast_dust(dust):
    too_fast = dust.select(lambda v: v.norm() >.05 | units.AUd,["velocity"])
    dust.remove_particles(too_fast)

if __name__ == "__main__":

    nstars = int(sys.argv[1])
    workers = int(sys.argv[2])
    method = sys.argv[3]
    print nstars
    seed = None

    n_dust = nstars

    convert_nbody = nbody_system.nbody_to_si(5.9736e24 | units.kg, 
                                             6371 | units.km)

    stars = core.Stars(n_dust)
    generate_dust(stars)
    remove_core_dust(stars)
    
    sun = stars[0]
    sun.mass = 1 | units.MSun
    sun.position = units.m(np.array((0.0,0.0,0.0)))
    sun.velocity = units.ms(np.array((0.0,0.0,0.0)))
    sun.radius = 0.05 | units.km

    earth = stars[1]
    earth.mass = units.kg(5.9736e24)
    earth.position = units.AU(np.array((8.418982185410142E-01,  5.355823303978186E-01,  2.327960005926782E-05)))
    vector0=np.matrix([8.418982185410142E-01,  5.355823303978186E-01,  2.327960005926782E-05])
    earth.velocity = units.AUd(np.array((-9.488931818313919E-03,  1.447515189957170E-02,  3.617712172296458E-07)))
    earth.radius = 0.05|units.km

    s = Viewer(stars,(1024, 780), convert_nbody)

    #stars.scale_to_standard()
    if method == 'octgrav':
        gravity = Octgrav()
    elif method == 'ph4':
        gravity = ph4()
    elif method == 'phigrape':					
        gravity = PhiGRAPE(convert_nbody)
    elif method == 'bhtree':
        gravity = BHTree(convert_nbody,number_of_workes = workers) 
    elif method == 'hermite':
        gravity = Hermite(convert_nbody,
                          number_of_workers = workers
                          #debugger = "xterm",
                          #redirection = "none"
                          )
    gravity.initialize_code()
    #gravity.parameters.epsilon_squared = 0.1 | nbody_system.length ** 2
    gravity.parameters.epsilon_squared = 0.01 | units.AU ** 2
    
    #stars.radius = 0.000001 | nbody_system.length
    gravity.particles.add_particles(stars)
    #gravity.stopping_conditions.pair_detection.enable()
    from_model_to_gravity = stars.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(stars)
    
    gravity.commit_particles()
    model_time = 0.0
    pairs = 0
    image_count = 0
    while s.go:
        flag = time.time()
        #s.animate()
        model_time +=1

        gravity.evolve_model(model_time|units.day)
        from_gravity_to_model.copy()
        
        rterr = 0
        while (time.time() - flag)<0.001:
            rterr +=1
        if rterr == 1: 
            print "realtime error"

        #if gravity.stopping_conditions.pair_detection.is_set():
        #    print gravity.stopping_conditions.pair_detection.particles(0).key
        #    print gravity.stopping_conditions.pair_detection.particles(1).key
        #    stars.remove_particles(gravity.stopping_conditions.pair_detection.particles(1))

        remove_core_dust(stars)
        remove_fast_dust(stars)
        from_model_to_gravity.copy()

        image_count +=1
        vector1 = np.matrix([earth.position[0].value_in(units.AU), 
                             earth.position[1].value_in(units.AU), 
                             earth.position[2].value_in(units.AU)])
        
        s.gamma = diffangle(vector0, vector1)/3.1415*180
        print s.gamma
        s.renderamuse(stars, "pic%03d.jpg" % image_count)
        s.handle_events()
        
    gravity.stop()
