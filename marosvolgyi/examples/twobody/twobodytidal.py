from amuse.community import *
import cPickle

import sys
sys.path.append('../../twobodymp')
import tidal_field
from interface import twobodympInterface
from interface import twobodymp

from amuse.units.units import *
from amuse.units import nbody_system as NU

import numpy
import time
from mpmath import *

PRECISION = 17

convert1 = NU.nbody_to_si(10**9 | MSun, 1000 | parsec)
convert2 = NU.nbody_to_si(1 | MSun, 100000 | AU)

convert_1_2=lambda x : convert2.to_nbody(convert1.to_si(x)) 

def total_energy_perturbed(nb, perturber):
  state=nb.get_state(0)
  time,err=nb.get_time()
  Ep,err=nb.get_potential_energy()
  Ek,err=nb.get_kinetic_energy()
  result,err=perturber.get_potential_at_point(0.,state['x'],state['y'],state['z'])  
  Es=result
  Etot=Ep+Es+Ek
  print "time = %.5f, Ep = %.6f, Ek = %.6f, Ep = %.6f, Etot = %.6f" % \
      (time,Ep,Ek,Es,Etot)

def tidal_field_from_file(fnaam):
    f=open(fnaam,'r')
    time=convert_1_2( (NU.time).new_quantity(cPickle.load(f)) )
    tides=convert_1_2( (NU.time**-2).new_quantity(cPickle.load(f)) )
    f.close() 
    tf=tidal_field.time_dependent_tidal_field(time.number,tides.number)
    return tf

def kick(instance, dtmodel):
    x0mp = "".zfill(120)
    y0mp = "".zfill(120)
    z0mp = "".zfill(120)
    vx0mp = "".zfill(120)
    vy0mp = "".zfill(120)
    vz0mp = "".zfill(120)
    x1mp = "".zfill(120)
    y1mp = "".zfill(120)
    z1mp = "".zfill(120)
    vx1mp = "".zfill(120)
    vy1mp = "".zfill(120)
    vz1mp = "".zfill(120)
    
    x0mp, y0mp, z0mp, res = instance.get_position_mp(0, x0mp,y0mp,z0mp, PRECISION)
    vx0mp, vy0mp, vz0mp, res = instance.get_velocity_mp(0, vx0mp,vy0mp,vz0mp, PRECISION)
    x1mp, y1mp, z1mp, res = instance.get_position_mp(1, x1mp,y1mp,z1mp, PRECISION)
    vx1mp, vy1mp, vz1mp, res = instance.get_velocity_mp(1, vx1mp,vy1mp,vz1mp, PRECISION)
    
    vx0 = mpf(vx0mp[0])#instance.get_velocity(0).vx_
    vy0 = mpf(vy0mp[0])#instance.get_velocity(0).vy_
    vz0 = mpf(vz0mp[0])#instance.get_velocity(0).vz_
    vx1 = mpf(vx1mp[0])#instance.get_velocity(0).vx_
    vy1 = mpf(vy1mp[0])#instance.get_velocity(0).vy_
    vz1 = mpf(vz1mp[0])#instance.get_velocity(0).vz_
    
    x0 = mpf(x0mp[0])
    y0 = mpf(y0mp[0])
    z0 = mpf(z0mp[0])
    x1 = mpf(x1mp[0])
    y1 = mpf(y1mp[0])
    z1 = mpf(z1mp[0])
    
    ax0, ay0, az0, error = get_gravity(1, x0, y0, z0)
    ax1, ay1, az1, error = get_gravity(1, x1, y1, z1)
    
    vx0 += dt_kick/2 * ax0
    vy0 += dt_kick/2 * ay0
    vz0 += dt_kick/2 * az0
    
    vx1 += dt_kick/2 * ax1
    vy1 += dt_kick/2 * ay1
    vz1 += dt_kick/2 * az1
    
    instance.set_velocity_mp(0, nstr(vx0,100), nstr(vy0,100), nstr(vz0,100))
    instance.set_velocity_mp(1, nstr(vx1,100), nstr(vy1,100), nstr(vz1,100))
    print x0,y0,z0,x1,y1,z1    
    
if __name__ == '__main__':
    perturber = tidal = tidal_field_from_file('../twobodymp/tideseq-1013801.dat')

    instance = twobodympInterface(redirection='none')
    instance.initialization(PRECISION)
    mp.dps = PRECISION

    mu=1.
    rinit=.2
    vinit=2.95
    
    energy=0.5*vinit**2-mu/rinit
    a=-mu/(2*energy)
    ecc=1-rinit/a
    rmax=a*(1+ecc)
    print a,ecc,rmax
    tperiod=2*numpy.pi/numpy.sqrt(mu)*a**(3./2)
    

    result = instance.new_particle([0.999, 0.001], [1,1],
                                   [0.0, rinit], [0.0, 0], [0.0, 0.0],
                                   [0.0, 0.0], [0.0, vinit], [0.0,0.0])

    
    dtmodel = 0.001*tperiod
    t = 0.0

    instance.commit_particles()
    instance.viewer(1)

    get_gravity = perturber.get_gravity_at_point

    for i in numpy.arange(0.000001, tperiod*10, dtmodel):
        instance.plot(150,50,50)
        tnow, err = instance.get_time()
        dt_kick = i-tnow
        print dtmodel, dt_kick

        perturber.evolve(tnow)
        kick(instance, dtmodel)
        instance.evolve_system(i)
        perturber.evolve(i)
        kick(instance, dtmodel)
        instance.plot(255,255,255)
  

        #total_energy_perturbed(instance, perturber)
        
    s=raw_input()

    instance.stop()
    
