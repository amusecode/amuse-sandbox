from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from interface import BonsaiInterface
#from interface import Bonsai

import os
import sys 
import numpy

from amuse.community.octgrav.interface import OctgravInterface, Octgrav

from amuse.units import nbody_system
from amuse.units import units
from amuse.support.codes import channel

from amuse.ext.plummer import *





from amuse import datamodel
class BonsaiInterfaceTests(TestWithMPI):
    
    def test0(self):
        print "Instantie aanmaken"
        instance = BonsaiInterface(redirection='none')
        print "aangemaakt"
        result,error = instance.echo_int(12)
        print "call instance done"
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
    def test1(self):
        plummer_size = 500
#        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        plummer =  MakePlummerModel(plummer_size)
        stars = plummer.result
#        stars.radius = range(1, plummer_size+1)|units.km
        mass=stars.mass.number
        radius=stars.radius.number
        x=stars.x.number
        y=stars.y.number
        z=stars.z.number
        vx=stars.vx.number
        vy=stars.vy.number
        vz=stars.vz.number


#        instance = self.new_instance_of_an_optional_code(Octgrav, convert_nbody)
        instance = BonsaiInterface(redirection='none') 
        instance.initialize_code()
        #ids,err = instance.new_particles(mass,radius,x,y,z,vx,vy,vz)
        for i in range(0,plummer_size):
            ids,err = instance.new_particle(mass[i],radius[i],x[i],y[i],z[i],vx[i],vy[i],vz[i])
            #print ids,err
        instance.commit_particles()
        instance.evolve_model(1)

        #energy_total_init = instance.potential_energy + instance.kinetic_energy
        #instance.evolve_model(100 | units.day)
        #energy_total_final = instance.potential_energy + instance.kinetic_energy

        #self.assertAlmostRelativeEqual(energy_total_init, energy_total_final, 2)
        mass,radius,x,y,z,vx,vy,vz,err=instance.get_state(ids)

        instance.stop()

