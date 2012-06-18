from amuse.couple import bridge

from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.smalln.interface import SmallN

from amuse.ic import plummer
from amuse.ic import gasplummer

from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units import nbody_system
from amuse.couple import multiples
from amuse import io

from optparse import OptionParser

import numpy
import time
import os.path
try:
    import pylab
except ImportError:
    pylab = None







class MultiplesRun(object):


    def setup_from_parameters(self,
        nstars = 10, 
        star_code = 'hermite',
        end_time = 10,
        seed = -1,
        snapshot_every = 1,
        star_radius = 0.1,
        star_smoothing_fraction = 0.0,
        base_name = 'test'
    ):
        
        if seed >= 0:
            numpy.random.seed(seed)
    
        self.snapshot_every = snapshot_every
        self.nstars = nstars
            
        self.time = 0 | nbody_system.time
        self.star_radius = star_radius | nbody_system.length
        self.star_epsilon = star_smoothing_fraction * self.star_radius
        self.end_time = end_time | nbody_system.time
        self.delta_time = snapshot_every | nbody_system.time
        
        dirname = ''
        for i in range(1000):
            dirname = '{0}{1:03d}'.format(base_name, i)
            if not os.path.exists(dirname):
                break
        
        self.runname = dirname
        print "self starting run:", self.runname
        
        os.makedirs(self.runname)
        
        self.particles = self.new_particles_cluster()
        print self.particles.velocity.lengths().mean()
        ke_mean = 0.5 * self.particles.mass.mean() * (self.particles.velocity.lengths().mean() ** 2)
        print  (nbody_system.G * (self.particles.mass.mean() ** 2) / 2 * ke_mean) / 2
        print ke_mean
        print nbody_system.energy
        
        self.code = getattr(self, 'new_star_code_{0}'.format(star_code))()
        self.code.particles.add_particles(self.new_particles_cluster())
        self.code.commit_particles()
        
        self.multiples = multiples.Multiples(self.code, self.new_smalln)
        
        self.snapshot()
        
        
    def new_smalln(self):
        result = SmallN()
        result.parameters.timestep_parameter = 0.1
        result.parameters.cm_index = 40000
        return result
        
    def new_particles_cluster(self):
        particles=plummer.new_plummer_model(self.nstars)
        particles.radius= self.star_radius
        particles.scale_to_standard()
        return particles

    def stop(self):
        pass
    
    def evolve_model(self):
            
        while self.time < self.end_time:
            self.time += self.delta_time
            
            self.multiples.evolve_model(self.time)
            
            self.snapshot()
                
    def snapshot(self):
        time = self.time.value_in(nbody_system.time)
        timestr = '{0:06.2f}'.format(time)
        
        subdir = os.path.join(self.runname, timestr)
        os.makedirs(subdir)
        
        self.current_snapshot = subdir
        
        #with open(os.path.join(subdir, 'parameters.pickle'), 'w') as f:
        #    pickle.dump(self, f)
        
        io.write_set_to_file(self.particles, os.path.join(subdir, 'particles.hdf5') , 'hdf5')
        for particle, tree in self.multiples.root_to_tree.iteritems():
            io.write_set_to_file(tree.get_tree_subset(), os.path.join(subdir, '{0}.hdf5'.format(particle.key)), 'hdf5')
            
        
    def __getstate__(self):
        parameters = {}
        
        return parameters
        
        
    def new_star_code_fi(self):
        result = Fi()
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = False
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.01 * self.delta_t
        return result
        
    def new_star_code_hermite(self):
        result = Hermite(number_of_workers=3)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        return result
        
    def new_star_code_phigrape(self):
        result = PhiGRAPE(mode="gpu")
        result.parameters.initialize_gpu_once = 1
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        return result
        
    def new_star_code_bhtree(self):
        result = BHTree()
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.01 * self.delta_t
        return result
        
    def new_star_code_octgrav(self):
        result = Octgrav()
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.01 * self.delta_t
        return result



        

def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n", "--nstar", 
        default = 10,
        dest="nstars",
        help="number of star particles",
        type="int"
    )
    result.add_option(
        "--star-code", 
        default = "hermite",
        dest="star_code",
        help="the code modelling the particles ('hermite', 'bhtree', 'octgrav', 'phigrape')",
        type="string"
    )
    
    result.add_option(
        "--star-smoothing-fraction", 
        default = 0.0,
        dest="star_smoothing_fraction",
        help="smoothing length of the stars as a fraction of the radii",
        type="float"
    )
    
    
    result.add_option(
        "--star-radius", 
        default = 0.01,
        dest="star_radius",
        help="radii of the stars, defaults to 0.01 nbody length units",
        type="float"
    )
    
    
    result.add_option(
        "-s", "--seed", 
        default = 0,
        dest="seed",
        help="random number seed (-1, no seed)",
        type="int"
    )
    result.add_option(
        "-t", "--end-time", 
        default = 10,
        dest="end_time",
        help="end time of the simulation (in nbody time, default 10)",
        type="float"
    )
    result.add_option(
        "--snapshot-every", 
        default = 1,
        dest="snapshot_every",
        help="delta timestep to snapshot (in nbody time, default 1)",
        type="float"
    )
    result.add_option(
        "--name", 
        default = "test",
        dest="base_name",
        help="name of the run, will be suffixed with a run number",
        type="string"
    )
    
    return result
    
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    code = MultiplesRun()
    code.setup_from_parameters(**options.__dict__)
    code.evolve_model()
    
