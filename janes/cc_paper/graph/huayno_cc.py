"""

Panels of particle distributions with different fractal dimension at 
different levels of the connected components subdivision. Plotted are 
f_dim=1.6 (left column), 2.3 (middle) and 3.0 (right) for the indicated 
level of the recursive connected component search. Note the difference 
between structured particle distribution and more smooth distributions: 
in the low f_dim case the CC algorithme tends to find different seperate 
sub systems, where as the smoother distribution fall apart in one large 
connected component and a number of single particles. The expectation is 
that the CC algorithm works best in the former case, while reducing to 
an SF split in the latter. Particles that are in singletons or binary 
components are removed in the sebsequent higher levels of the recursion.

"""
import numpy

from amuse.ic.fractalcluster import new_fractal_cluster_model
from amuse.ic.plummer import new_plummer_model

from amuse.units import nbody_system
from amuse.units import units,constants
from amuse.units.quantities import zero
from matplotlib import pyplot

from amuse.ext.basicgraph import Graph
from amuse.ext.basicgraph import ConnectedComponents

from amuse.datamodel import Particles

colors=["r","g","b","c","k","y"]
markers=["o","x","d","*","+"]

def timestep(ipart,jpart, dt_param=0.05, rarvratio=1.,_G=constants.G):
    eps=1.e-20 | zero.unit
    
    dx=ipart.x-jpart.x  
    dy=ipart.y-jpart.y
    dz=ipart.z-jpart.z
    dr2=dx**2+dy**2+dz**2+eps
#    if dr2>0:
    dr=dr2**0.5
    dr3=dr*dr2
    dvx=ipart.vx-jpart.vx  
    dvy=ipart.vy-jpart.vy  
    dvz=ipart.vz-jpart.vz
    vdotdr2=(dvx*dx+dvy*dy+dvz*dz)/dr2
    dv2=dvx**2+dvy**2+dvz**2+eps
    mu=_G*(ipart.mass+jpart.mass)
    tau=rarvratio*dt_param/2.**0.5*(dr3/mu)**0.5
    dtau=3*tau*vdotdr2/2.
    dtau=numpy.minimum(dtau,numpy.ones_like(dtau))
    tau1=tau/(1-dtau/2)          
#      if dv2>0:
    tau=dt_param*dr/dv2**0.5
    dtau=tau*vdotdr2*(1+mu/(dv2*dr))  
    dtau=numpy.minimum(dtau,numpy.ones_like(dtau))
    tau2=tau/(1-dtau/2)          
    return tau1.unit(numpy.minimum(tau1.value_in(tau1.unit),tau2.value_in(tau1.unit)))

def connected_components(parts, treshold=None, distfunc=None):
  if distfunc is None:
    def distfunc(p,q):
      return ((p.x-q.x)**2+(p.y-q.y)**2+(p.z-q.z)**2)**0.5

  print "making graph"
  graph=Graph()
  for p in parts:
     graph.add_node(p)
     d=distfunc(p,parts)
     edges= [ (d[i],p,q) for i,q in enumerate(parts) if p!=q] 
     edges=filter( lambda x: x[0] < treshold,edges)
     for e in edges:
       graph.add_edge(e[1],e[2],e[0]) 
  print "done"
  cc=ConnectedComponents(graph)
  print "number of edges:", len(graph.all_edges())
  print "number of CC:",len(cc)
  return graph,cc

def testcc():  
  N=200
  
  parts = new_fractal_cluster_model(N=N,fractal_dimension=2.0,random_seed=1234532)
#  parts = new_plummer_model(N)
    
  f=pyplot.figure(figsize=(8,8))
  pyplot.plot(parts.x.number,parts.y.number,'r.')
  pyplot.xlim(-1,1)
  pyplot.ylim(-1,1)
  pyplot.savefig("fractal.png")
      
  def distfunc(p,q):    
    return timestep(p,q,_G=nbody_system.G)
        
  graph,cc=connected_components(parts,treshold=1./512 | nbody_system.time,distfunc=distfunc)
  
  f=pyplot.figure(figsize=(8,8))  
  for i,parts in enumerate(cc):
    for p in parts:    
      pyplot.plot(p.x.number,p.y.number,colors[i%len(colors)]+'.')
  alledges=graph.all_edges()
  for e in alledges:
    pyplot.plot([e[2].x.number,e[1].x.number],[e[2].y.number,e[1].y.number],'b')
  pyplot.xlim(-1,1)
  pyplot.ylim(-1,1)
  pyplot.savefig("cc.png")

def testhuayno(N=100,fd=1.6,seed=1234567):
  
  parts = new_fractal_cluster_model(N=N,fractal_dimension=fd,random_seed=seed)

#  parts.z*=0
#  parts.vz*=0

  def distfunc(p,q):    
    return timestep(p,q,_G=nbody_system.G)

  f=pyplot.figure(figsize=(8,8))
  pyplot.plot(parts.x.number,parts.y.number,'r.')
  pyplot.xlim(-1,1)
  pyplot.ylim(-1,1)
  pyplot.savefig("fractal-%3.1f.png"%fd)

  i=0
  cc=[]
  newparts=parts
  while len(cc)<len(parts):
    tcurrent=1./2**i |  nbody_system.time
#    lcurrent=1./2**i |  nbody_system.length
    parts=newparts
    graph,cc=connected_components(parts,treshold=tcurrent,distfunc=distfunc)
#    graph,cc=connected_components(parts,treshold=lcurrent)

    print i,tcurrent,len(cc),len(parts)

    if len(cc) > 1:
      f=pyplot.figure(figsize=(8,8))  
      alledges=graph.all_edges()
      for e in alledges:
        pyplot.plot([e[2].x.number,e[1].x.number],[e[2].y.number,e[1].y.number],'grey',linewidth=0.5)
      for ip,iparts in enumerate(cc):
        for p in iparts:    
          pyplot.plot(p.x.number,p.y.number,'k.',markersize=8.,mew=.5)
      pyplot.xlim(-1.2,1.2)
      pyplot.ylim(-1.2,1.2)

      single=0
      for s in cc:
        if len(s)==1:
          single+=1

      binary=0
      for s in cc:
        if len(s)==2:
          binary+=1

      pyplot.text(0.6,-.85,"level: %2i"%i,fontsize=16)
      pyplot.text(0.6,-.95,"#cc: %i"%len(cc),fontsize=16)
      pyplot.text(0.6,-1.05,"#s,b: %i,%i"%(single,binary),fontsize=16)
      pyplot.text(0.6,-1.15,"#parts: %i"%len(parts),fontsize=16)

      pyplot.savefig("cc-%3.1f-%3.3i.png"%(fd,i))
      
      newparts=Particles()
      for s in cc:
        if len(s)>2:
          for p in s:
            newparts.add_particle(p)

    i+=1
  print len(cc),len(parts)

if __name__=="__main__":
#  import cProfile
#  cProfile.run("timecc()","prof")
  testhuayno(1000,1.6,seed=1234561)
  testhuayno(1000,2.3,seed=1234563)
  testhuayno(1000,3.0,seed=1234561)
