from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import pyopencl as cl
import numpy as np
from amuse.ext.plummer import MakePlummerModel
import sys
import ctypes
import time

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
nbodies = 4096
frame=0
t0=time.time()
t1=t0

parts = MakePlummerModel(nbodies).result

#x_r = np.random.rand(nbodies,4).astype('float32')
#v_r = 0.0 * np.ones((nbodies,4)).astype('float32')

x_r = np.vstack([parts.x.number,parts.y.number,parts.z.number, parts.mass.number]).transpose()
x_r = np.array(x_r,dtype='float32')
x_r = x_r.flatten()
v_r = np.vstack([parts.vx.number,parts.vy.number,parts.vz.number, np.zeros(nbodies)]).transpose()
v_r = np.array(v_r,dtype='float32')
v_r = v_r.flatten()

ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

mf = cl.mem_flags

x_buf = cl.Buffer(ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=x_r)
v_buf = cl.Buffer(ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=v_r)
x_re = cl.Buffer(ctx, mf.READ_WRITE , x_r.nbytes)
v_re = cl.Buffer(ctx, mf.READ_WRITE , v_r.nbytes)

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
    float4 a= (float4) (0.0f,0.0f,0.0f,0.0f);
    float invdist;
    float invdist3;
    const float4 dt1 = (float4)(dt,dt,dt,0.0f);


    for (i = 0; i< nbodies; i++) {
      dist = r[i] - r[gid];
      invdist = rsqrt(dist.x*dist.x + dist.y*dist.y + dist.z*dist.z + eps2);
      invdist3 = r[i].w*invdist*invdist*invdist;
      float4 F= (float4) (invdist3,invdist3,invdist3,0.0f);
      a += invdist3*dist;
   }
   vr[gid] = v[gid] + dt1 * a ;
   rr[gid] = r[gid] + dt1 * vr[gid];

}""")

try:
    prg.build()
except:
    print("Error:")
    print(prg.get_build_info(ctx.devices[0], cl.program_build_info.LOG))
    raise

t0=time.time()
t1=t0

def display():

    global x_buf
    global v_buf
    global x_re
    global v_re
    global zoom
    global x_r
    global v_r
    global  frame
    global t0
    global t1
    
    prg.demo(queue, (nbodies,), None, x_buf, v_buf, x_re, v_re, 
             np.array([0.001],dtype='float32'), 
             np.array([0.001],dtype='float32'), 
             np.array([nbodies],dtype='int32'))
    prg.demo(queue, (nbodies,), None, x_re, v_re, x_buf, v_buf, 
             np.array([0.001],dtype='float32'), 
             np.array([0.001],dtype='float32'), 
             np.array([nbodies],dtype='int32'))
    
    cl.enqueue_read_buffer(queue, x_buf, x_r).wait()
#    cl.enqueue_read_buffer(queue, v_buf,  v_r).wait()

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
    
    glEnableClientState(GL_VERTEX_ARRAY)
#    glBegin(GL_POINTS)
#    print x_r.shape
    p=x_r.reshape((nbodies,4))
#    print p 
    x=np.array(p[0:nbodies,0:3],dtype='float32')

#    print x

    glVertexPointerf(x)
#    for xx in x:
#      glVertex3f(xx[0],xx[1],xx[2])
    

#    glEnd()
    glDrawArrays(GL_POINTS,0,nbodies)

    glColor4f(0.0, 0.0, 1.0, 1.0)
    glutWireCube(1.0)
    glutWireCube(2.0)
    glutWireCube(4.0)
    glutWireCube(8.0)

    glPopMatrix()
    glutSwapBuffers()
    frame+=1
    
    if frame%10==0:
      t=time.time()
      print 'fps:',frame/(t-t0),'  gflops:', 100*2.*nbodies**2*20/(t-t0)/1.e9
      t0=t
      frame=0

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
