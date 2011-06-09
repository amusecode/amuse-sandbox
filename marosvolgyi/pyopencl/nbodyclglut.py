from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import pyopencl as cl
import numpy as np
from amuse.ext.plummer import MakePlummerModel
import sys

drag = False
trail = False
box = False
alpha = 0.0
beta = 0.0
gamma = 0.0
omega = 0.0
zoom = -12
xold = 0
yold = 0
nbodies = 400

parts = MakePlummerModel(nbodies).result

#x_r = np.random.rand(nbodies,4).astype('float32')
#v_r = 0.0 * np.ones((nbodies,4)).astype('float32')

x_r = np.array([parts.x.number,parts.y.number,parts.z.number],
                dtype='float32').transpose()
print x_r
v_r = np.array([parts.vx.number,parts.vy.number,parts.vz.number],
                dtype='float32').transpose()
x_r=np.append(x_r, 1.0/nbodies * np.ones((nbodies,1)).astype('float32'),1)
v_r=np.append(v_r, np.zeros((nbodies,1)).astype('float32'),1)
print x_r;#exit()
ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

mf = cl.mem_flags

x_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=x_r)
v_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=v_r)
x_re = cl.Buffer(ctx, mf.WRITE_ONLY , x_r.nbytes)
v_re = cl.Buffer(ctx, mf.WRITE_ONLY , v_r.nbytes)

prg = cl.Program(ctx,
"""
__kernel void demo(__global float4 *r, 
                   __global float4 *v,
                   __global float4 *rr, 
                   __global float4 *vr,
                   float dt,
                   float eps2,
                   int nbodies)
{
    int i;
    int gid = get_global_id(0);
    float4 dist;
    float4 F;
    F.x = 0.0;F.y=0.0;F.z=0.0;F.w=0.0;
    float invdist2;
    float invdist3;

    for (i = 0; i< nbodies; ++i) {
      dist = r[gid] -r[i];
      invdist2 = rsqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z + eps2);
      invdist3 = invdist2*invdist2*invdist2;
      F -= .005 * dist * invdist3;
   }

   v[gid] += F * dt;
   r[gid] += v[gid] * dt;

   rr[gid] = r[gid];
   vr[gid] = v[gid];

}""")

try:
    prg.build()
except:
    print("Error:")
    print(prg.get_build_info(ctx.devices[0], cl.program_build_info.LOG))
    raise

def display():

    global x_buf
    global v_buf
    global x_re
    global v_re
    global zoom

    prg.demo(queue, (nbodies,4), None, x_buf, v_buf, x_re, v_re, 
             np.array([0.001],dtype='float32'), 
             np.array([0.01],dtype='float32'), 
             np.array([nbodies],dtype='int32'))
    cl.enqueue_read_buffer(queue, x_re,  x_r).wait()
    cl.enqueue_read_buffer(queue, v_re,  v_r).wait()

    x_buf = x_re
    v_buf = v_re

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    
    glPushMatrix()
    glEnable(GL_DEPTH_TEST)

    glLoadIdentity()
    gluPerspective( 20.0, 1.0, 1.0, 300.0)
    glClearColor(0.0 ,0.0, 0.0, 0.0)
    glClear(GL_COLOR_BUFFER_BIT)
    glClear(GL_DEPTH_BUFFER_BIT)

    glPointSize(2.0)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE)
    glEnable(GL_BLEND)
    glDepthMask(GL_FALSE)

    glTranslatef(0.0, 0.0, zoom)    

    glRotatef(beta, 1.0, 0.0, 0.0)
    glRotatef(alpha, 0.0, 1.0, 0.0)
    glRotatef(gamma, 0.0, 0.0, 1.0)
 
    glTranslate(0.0,0.0,+12)
	
	
    glColor4f(1.0, 1.0, 0.0, 0.5)
    glBegin(GL_POINTS)
    for i, r in enumerate(x_r):
        glVertex3f(r[0],r[1],r[2])
    glEnd()

    glColor4f(0.0, 0.0, 1.0, 1.0)
    glutWireCube(1.0)
    glutWireCube(2.0)
    glutWireCube(4.0)
    glutWireCube(8.0)

    glPopMatrix()
    glutSwapBuffers()


def init():
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);

    glEnable(GL_DEPTH_TEST)
    glMatrixMode(GL_PROJECTION)
    gluPerspective( 40.0, 1.0, 1.0, 300.0)
    glMatrixMode(GL_MODELVIEW)
    gluLookAt(0.0, 0.0, 10.0,  
              0.0, 0.0, 0.0,  
 	      0.0, 1.0, 0.)
	
#eye is at (0,0,5) */		        
#center is at (0,0,0) */		      
#up is in positive Y direction */      

def reshape(w, h):
   glViewport(0, 0, w, h)
   glMatrixMode (GL_PROJECTION)
   glLoadIdentity()
   gluPerspective( 40.0, w/h, 1.0, 100.0)
   gluLookAt(0.0, 0.0, 10.0,  
             0.0, 0.0, 0.0,  
             0.0, 1.0, 0.)

def keyboard(key, x, y):
    global zoom
    if key == chr(27):
        sys.exit(0)
    if (key == 'z'): zoom += 10
    if (key == 'Z'): zoom -= 10;
    print zoom

def mouseMoved(x, y):
    global xold 
    global yold  
    global alpha
    global beta
    dx = x - xold;
    dy = y - yold;
    if (dx > 0): alpha += .5
    if (dx < 0): alpha -= .5
    if (dy > 0): beta += .5
    if (dy < 0): beta -= .5
    xold = x;
    yold = y;

def mouseWheel(button, dire, x, y):
    global zoom
    print dire
    if (dire>0): zoom +=10
    if (dire<0): zoom -=10

if __name__ == '__main__':
    glutInit(sys.argv)
    glutInitDisplayMode (GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
    glutInitWindowSize (500, 500)
    glutCreateWindow('scene')
    init()
    glutReshapeFunc(reshape)
    glutKeyboardFunc(keyboard)
    glutMotionFunc(mouseMoved)  
    glutMouseWheelFunc(mouseWheel)
    glutIdleFunc(display)
    glutDisplayFunc(display)
    glutMainLoop()
