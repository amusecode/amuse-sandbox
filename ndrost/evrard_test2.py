"""
initial conditions for the SPH evrard collapse test
"""
from amuse.units import units
from amuse.units import nbody_system

from math import *
import numpy
from amuse.units import generic_unit_system as generic_system 
from amuse.units import nbody_system
from amuse import datamodel
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesWithUnitsConverted
class uniform_random_unit_cube(object):
  def __init__(self,targetN):
    self.targetN=targetN
    self.par=long(targetN)
  def make_xyz(self):
    x=numpy.random.uniform(-1.,1.,self.par)
    y=numpy.random.uniform(-1.,1.,self.par)
    z=numpy.random.uniform(-1.,1.,self.par)
    return x,y,z

class regular_grid_unit_cube(object):
  def __init__(self,targetN):
    self.targetN=targetN
    self.par=long(float(targetN)**(1./3.)+1.5) 
  def make_xyz(self):
    nf=self.par
    x,y,z=numpy.mgrid[-1.:1.:nf*1j,-1.:1.:nf*1j,-1.:1.:nf*1j] 
    x=x.flatten()
    y=y.flatten()
    z=z.flatten()
    return x,y,z

class body_centered_grid_unit_cube(object):
  def __init__(self,targetN):
    self.targetN=targetN
    self.par=long(float(targetN/2.)**(1./3.)+1.5) 

  def make_xyz(self):
    nf=self.par
    x1,y1,z1=numpy.mgrid[-1.:1.:nf*1j,-1.:1.:nf*1j,-1.:1.:nf*1j] 
    x2,y2,z2=numpy.mgrid[-1.+1./2/nf:1.-1./2/nf:(nf-1)*1j, \
      -1.+1./2/nf:1.-1./2/nf:(nf-1)*1j,-1.+1./2/nf:1.-1./2/nf:(nf-1)*1j]                        
    x=numpy.concatenate( (x1.flatten(),x2.flatten()) )
    y=numpy.concatenate( (y1.flatten(),y2.flatten()) )
    z=numpy.concatenate( (z1.flatten(),z2.flatten()) )
    return x,y,z

class uniform_unit_sphere(object):
  def __init__(self,targetN, base_grid=None):
    cube_sphere_ratio=4/3.*numpy.pi*0.5**3
    self.targetN=targetN/cube_sphere_ratio
    if base_grid is None: 
      self.base_grid=uniform_random_unit_cube(self.targetN)
    else:
      self.base_grid=base_grid(self.targetN)
 
  def cutout_sphere(self,x,y,z):
    r=x**2+y**2+z**2
    selection=r < numpy.ones_like(r)        
    x=x.compress(selection)
    y=y.compress(selection)
    z=z.compress(selection)
    return x,y,z

  def make_xyz(self):
    return self.cutout_sphere(*self.base_grid.make_xyz())
        
class MakeEvrardTest(object):
    def __init__(self, targetN, convert_nbody = None, base_grid=None, size=1.,
                   mass=1.,internal_energy=0.05,seed=345672):
        numpy.random.seed(seed)
        self.targetN = targetN
        self.convert_nbody = convert_nbody
        self.size=size
        self.mass=mass
        self.internal_energy=internal_energy
        self.base_sphere=uniform_unit_sphere(targetN,base_grid)   
           
    def new_model(self):
        x,y,z=self.base_sphere.make_xyz()
        self.actualN=len(x)
        r=numpy.sqrt(x**2+y**2+z**2)
        rtarget=self.size*r**1.5
        mass=numpy.ones_like(x)/self.actualN
        internal_energy=numpy.ones_like(x)*self.internal_energy
        r=r.clip(1.e-8,2*self.size)
        x=rtarget*x/r
        y=rtarget*y/r
        z=rtarget*z/r
        vx=numpy.zeros_like(x)
        vy=numpy.zeros_like(x)
        vz=numpy.zeros_like(x)
        return (mass,x,y,z,vx,vy,vz,internal_energy)
    
    @property
    def result(self):
        mass,x,y,z,vx,vy,vz,u = self.new_model()
        result = datamodel.Particles(self.actualN)
        result.mass = nbody_system.mass.new_quantity(mass)
        result.x = nbody_system.length.new_quantity(x)
        result.y = nbody_system.length.new_quantity(y)
        result.z = nbody_system.length.new_quantity(z)
        result.vx = nbody_system.speed.new_quantity(vx)
        result.vy = nbody_system.speed.new_quantity(vy)
        result.vz = nbody_system.speed.new_quantity(vz)
        result.u = (nbody_system.speed**2).new_quantity(u)

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_generic())
            result = result.copy_to_memory()
            
        return result


if __name__=="__main__":
  convert_nbody = nbody_system.nbody_to_si(100. | units.MSun, 1.0 | units.parsec)
  sphere=MakeEvrardTest(10000,convert_nbody)
  parts=sphere.result
  print parts[0].internal_energy**0.5
  print len(parts)
  
