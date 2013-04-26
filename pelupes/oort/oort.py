import numpy

from amuse.units import units,constants

from amuse.community.mercury.interface import Mercury

from amuse.datamodel import Particles

from matplotlib import pyplot

ss=Particles(4)

ss[0].mass=1.| units.MSun
ss[0].radius=1.| units.RSun

mu=ss[0].mass*constants.G
rmax=[10000.,20000.,50000.] | units.AU 
eps=0.5
rmin=(1-eps)/(1+eps)*rmax
a=(rmin+rmax)/2
v0=mu**0.5*(2./rmax-1/a)**0.5
T=2*numpy.pi*(a**3/mu)**0.5

ss[1:4].radius=1.| units.km
ss[1:4].x=rmax
ss[1:4].y=0 | units.AU
ss[1:4].z=0 | units.AU
ss[1:4].vx=0. | units.kms
ss[1:4].vy=v0.in_(units.kms)
ss[1:4].vz=0. | units.kms

grav=Mercury()

grav.parameters.timestep=T[0]/400

grav.particles.add_particles(ss)

tnow=grav.model_time
tend=T[-1]

data=dict()
data['x']=[]
data['y']=[]

while tnow<tend:
  tnow+=T[0]/400.*8
  grav.evolve_model(tnow)
  data['x'].append(grav.particles.x[1:4].value_in(units.AU))
  data['y'].append(grav.particles.y[1:4].value_in(units.AU))

xx=numpy.array(data['x'])
yy=numpy.array(data['y'])

pyplot.plot(xx[:,0],yy[:,0])
pyplot.plot(xx[:,1],yy[:,1])
pyplot.plot(xx[:,2],yy[:,2])
pyplot.xlim(-50000,50000)
pyplot.ylim(-50000,50000)
pyplot.show()
