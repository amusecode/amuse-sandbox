from amuse.couple import bridge

from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.phiGRAPE.interface import PhiGRAPE

from amuse.ic import plummer
from amuse.ic import gasplummer

from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units import nbody_system

from optparse import OptionParser

import numpy
import time
try:
	import pylab
except ImportError:
	pylab = None







class MultiplesRun(object):


    def setup_from_parameters(self,
        nstars = 10, 
        total_mass = -1,
        end_time = 10,
        seed = -1,
        snapshot_every = 1,
        star_radius = 0.1,
        star_smoothing_fraction = 0.001,
        basename = 'test'
    ):
        
        if seed >= 0:
            numpy.random.seed(seed)
    
        self.snapshot_every = snapshot_every
        self.nstars = nstars
        
        if total_mass <= 0.0:
            total_mass = 1
            
        self.time = 0 | nbody_system.time
        self.total_mass = total_mass | nbody_system.mass
        self.star_radius = star_radius | nbody_system.length
        self.star_epsilon = star_smoothing_fraction * self.star_radius
        self.end_time = end_time | nbody_system.time
        self.delta_t = snapshot_every | nbody_system.time
        
        dirname = ''
        for i in range(100):
            dirname = '{0}{1:03i}'.format(basename, i)
            if not os.path.exists(dirname):
                break
        
        self.runname = dirname
        print "self starting run:", self.runname
        
        os.mkdir(self.runname)
        
        self.particles = self.new_particles_cluster()
        self.snapshot()
        
    def new_particles_cluster(self):
        particles=plummer.new_plummer_sphere(self.nstars,convert_nbody=self.converter)
        particles.radius= self.star_epsilon
        return particles

    def stop(self):
        pass
    
    def evolve_model(self):
            
        while self.time < self.end_time
            self.time += self.delta_time
            
            self.code.evolve_model(self.time)
            
            print self.code.time
            
            self.snapshot()
                
    def snapshot(self):
        time = self.time.value_in(nbody_system.time)
        timestr = '{0:6.2f}'.format(time)
        
        subdir = os.path.join(self.runname, timestr)
        os.path.mkdir(subdir)
        
        self.snapshot = subdir
        
        with open(os.path.join(subdir, 'particles.hdf5')) as f:
            pickle.dump(self, f)
        
        write_set_to_file(self.particles, os.path.join(subdir, 'particles.hdf5') , 'hdf5')
        for particle, tree in self.code.root_to_tree.iteritems():
            write_set_to_file(tree, os.path.join(subdir, '{0}.hdf5'.format(particle.key)), 'hdf5')
            
        
        
        
        


class BridgeStarAndGasPlummerCode(AbstractStarAndGasPlummerCode):


    def __init__(self,
        nstars = 10, 
        ngas = -1, 
        endtime = 10,
        total_mass = 1000,
        gas_fraction = 0.9,
        rscale = 1.0,
        star_code = 'hermite',
        gas_code = 'field', 
        star_smoothing_fraction = 0.001,
        gas_smoothing_fraction = 0.05,
        seed = -1,
        ntimesteps = 10,
        interaction_timestep = 0.01,
        must_do_plot = True,
        gas_to_star_interaction_code = 'none',
        star_to_gas_interaction_code = 'none',
        **ignored_options
    ):
        
        AbstractStarAndGasPlummerCode.__init__(
            self,
            nstars, 
            ngas, 
            endtime,
            total_mass,
            gas_fraction,
            rscale,
            star_smoothing_fraction,
            gas_smoothing_fraction,
            seed,
            ntimesteps,
            must_do_plot
        )

        
        self.interaction_timestep = self.converter.to_si(interaction_timestep| nbody_system.time)
        
        self.create_codes(
            gas_code,
            star_code,
            gas_to_star_interaction_code,
            star_to_gas_interaction_code,
        )
        
        self.create_bridge()
        
        self.code = self.bridge_system
        
        time = 0
        sum_energy = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        energy = self.converter.to_nbody(sum_energy).value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())

        print "Time          :", time
        print "Energy        :", energy
        print "Virial radius :", coreradius
        
        self.evolve_model()
        
        if must_do_plot:
            pylab.show()
            pylab.savefig(
                "{0}-{1}-{2}-{3}.png".format(
                    star_code,
                    gas_code,
                    nstars,
                    ngas
                )
            )
        
        
        time = self.converter.to_nbody(self.code.time).value_in(nbody_system.time)
        sum_energy = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        energy = self.converter.to_nbody(sum_energy).value_in(nbody_system.energy)
        coreradius = self.star_code.particles.virial_radius().value_in(self.rscale.to_unit())
        
        print "Time          :", time
        print "Energy        :", energy
        print "Virial radius :", coreradius
        
        self.stop()
        
        if must_do_plot:
            raw_input('Press enter...') 
        
   
        
    def create_codes(self, gas_code, star_code, gas_to_star_interaction_code, star_to_gas_interaction_code):
        self.star_code = getattr(self,'new_star_code_'+star_code)()
        self.gas_code = getattr(self, 'new_gas_code_'+gas_code)()
        
        self.gas_to_star_codes = getattr(self, 'new_gas_to_star_interaction_codes_'+gas_to_star_interaction_code)(self.gas_code)
        self.star_to_gas_codes = getattr(self, 'new_star_to_gas_interaction_codes_'+star_to_gas_interaction_code)(self.star_code)
        
    def create_bridge(self):
        bridge_code1 = bridge.GravityCodeInField(
            self.gas_code, self.star_to_gas_codes
        )
        bridge_code2 = bridge.GravityCodeInField(
            self.star_code, self.gas_to_star_codes
        )
        
        self.bridge_system = bridge.Bridge(
            timestep = self.interaction_timestep,
            use_threading = False
        )
        self.bridge_system.add_code(bridge_code2)
        self.bridge_system.add_code(bridge_code1)
    
    def stop(self):
        self.star_code.stop()
        self.gas_code.stop()
    
    def new_gas_to_star_interaction_codes_self(self, gas_code):
        return [gas_code]
        
    def new_star_to_gas_interaction_codes_self(self, star_code):
        return [star_code]
        
    def new_gas_to_star_interaction_codes_none(self, gas_code):
        return []
        
    def new_star_to_gas_interaction_codes_none(self, gas_code):
        return []
        
    def new_gas_to_star_interaction_codes_octgrav(self, gas_code):
        def new_octgrav():
            result = Octgrav(self.converter)
            result.parameters.epsilon_squared = self.gas_epsilon ** 2
            return result
            
        return [bridge.CalculateFieldForCodes(new_octgrav, [gas_code])]
    
    
    def new_gas_to_star_interaction_codes_bhtree(self, gas_code):
        def new_bhtree():
            result = BHTree(self.converter)
            result.parameters.epsilon_squared = self.star_epsilon ** 2
            return result
            
        return [bridge.CalculateFieldForCodes(new_bhtree, [gas_code])]\
    
    def new_gas_code_fi(self):
        result = Fi(self.converter)
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = True
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.125 * self.interaction_timestep
        #result.parameters.adaptive_smoothing_flag = True
        #result.parameters.epsilon_squared = self.gas_epsilon ** 2
        #result.parameters.eps_is_h_flag = False
        result.parameters.integrate_entropy_flag = False
        
        #result.parameters.self_gravity_flag = False
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_fi(self):
        result = Fi(self.converter)
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = False
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_gas_code_gadget(self):
        result = Gadget2(self.converter)
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        
    def new_gas_code_field(self):
        result = GasPlummerModelExternalField(
            radius = self.rscale,
            total_mass = self.gas_mass
        )
        return result
        
    def new_gas_code_hermite(self):
        result = Hermite(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster_as_gas())
        result.commit_particles()
        return result
        
    def new_star_code_hermite(self):
        result = Hermite(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_phigrape(self):
        result = PhiGRAPE(self.converter, mode="gpu")
        result.parameters.initialize_gpu_once = 1
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_bhtree(self):
        result = BHTree(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_star_code_octgrav(self):
        result = Octgrav(self.converter)
        result.parameters.epsilon_squared = self.star_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster())
        result.commit_particles()
        return result
        
    def new_gas_code_bhtree(self):
        result = BHTree(self.converter)
        result.parameters.epsilon_squared = self.gas_epsilon ** 2
        result.parameters.timestep = 0.125 * self.interaction_timestep
        result.particles.add_particles(self.new_particles_cluster_as_gas())
        result.commit_particles()
        return result



class AllInOneStarAndGasPlummerCode(AbstractStarAndGasPlummerCode):


    def __init__(self,
        nstars = 10, 
        ngas = -1, 
        endtime = 10,
        total_mass = 1000,
        gas_fraction = 0.9,
        rscale = 1.0,
        sph_code = 'fi',
        star_smoothing_fraction = 0.001,
        gas_smoothing_fraction = 0.05,
        seed = -1,
        ntimesteps = 10,
        must_do_plot = True,
        interaction_timestep = 0.01,
        **ignored_options
    ):
        
        AbstractStarAndGasPlummerCode.__init__(
            self,
            nstars, 
            ngas, 
            endtime,
            total_mass,
            gas_fraction,
            rscale,
            star_smoothing_fraction,
            gas_smoothing_fraction,
            seed,
            ntimesteps,
            must_do_plot
        )

        
        self.interaction_timestep = self.converter.to_si(interaction_timestep| nbody_system.time)
        
        self.create_code(sph_code)
        
        sum_energy = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        energy = self.converter.to_nbody(sum_energy).value_in(nbody_system.energy)
        coreradius = self.code.dm_particles.virial_radius().value_in(self.rscale.to_unit())
        
        print "Time:", 0
        print "Energy:", energy
        print "Virial radius:", coreradius
        
        self.evolve_model()
        
        if must_do_plot:
            pylab.show()
            pylab.savefig(
                "{0}-{1}-{2}.png".format(
                    sph_code,
                    nstars,
                    ngas
                )
            )
        
        time = self.converter.to_nbody(self.code.model_time).value_in(nbody_system.time)
        sum_energy = self.code.kinetic_energy + self.code.potential_energy + self.code.thermal_energy
        energy = self.converter.to_nbody(sum_energy).value_in(nbody_system.energy)
        coreradius = self.code.dm_particles.virial_radius().value_in(self.rscale.to_unit())
        
        print "Time:", time
        print "Energy:", energy
        print "Virial radius:", coreradius
        
        self.stop()
        
        if must_do_plot:
            raw_input('Press enter...') 
        
    def evolve_model(self):
        
        if self.must_do_plot:
            self.update_plot(time = 0 * self.delta_t, code = self.code)
            
        for time in self.delta_t * range(1,self.ntimesteps+1):
            self.code.evolve_model(time)
            print self.converter.to_nbody(self.code.model_time)
            if self.must_do_plot:
                self.update_plot(time = self.code.time, code = self.code)
        
        
        
    def create_code(self, name):
        self.code = getattr(self, 'new_sph_code_'+name)()
        
    def stop(self):
        self.code.stop()
    
    def new_sph_code_fi(self):
        result = Fi(self.converter)
        result.parameters.self_gravity_flag = True
        result.parameters.use_hydro_flag = True
        result.parameters.radiation_flag = False
        result.parameters.periodic_box_size = 500 | units.parsec
        result.parameters.timestep = 0.125 * self.interaction_timestep
        
        #result.parameters.adaptive_smoothing_flag = True
        #result.parameters.epsilon_squared = self.gas_epsilon ** 2
        #result.parameters.eps_is_h_flag = False
        
        result.parameters.integrate_entropy_flag = False
        result.dm_particles.add_particles(self.new_particles_cluster())
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
        return result
        
    def new_sph_code_gadget(self):
        result = Gadget2(self.converter)
        result.dm_particles.add_particles(self.new_particles_cluster())
        result.gas_particles.add_particles(self.new_gas_cluster())
        result.commit_particles()
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
        "-g", "--ngas", 
        default = -1,
        dest="ngas",
        help="number of gas particles (if -1, 10 times the number of stars)",
        type="int"
    )
    result.add_option(
        "--gas-code", 
        default = "field",
        dest="gas_code",
        help="the code modelling the gas ('fi', 'gadget', 'field')",
        type="string"
    )
    result.add_option(
        "--star-code", 
        default = "hermite",
        dest="star_code",
        help="the code modelling the particles ('hermite', 'bhtree', 'octgrav', 'phigrape')",
        type="string"
    )
    result.add_option(
        "--sph-code", 
        default = "fi",
        dest="sph_code",
        help="the code modelling the particles and the gas simultaniously",
        type="string"
    )
    result.add_option(
        "--gas-star-code", 
        default = "self",
        dest="gas_to_star_interaction_code",
        help="the code calculating the gravity field of the gas code for the star code (default is self, gas code will calculate field for star code)",
        type="string"
    )
    result.add_option(
        "--star-gas-code", 
        default = "self",
        dest="star_to_gas_interaction_code",
        help="the code calculating the gravity field of the star code for the gas code (default is self, star code will calculate field for gas code)",
        type="string"
    )
    result.add_option(
        "-m", "--total-mass", 
        default = 1000.0,
        dest="total_mass",
        help="the total mass in solar masses",
        type="float"
    )
    result.add_option(
        "--gas-fraction", 
        default = 0.9,
        dest="gas_fraction",
        help="the gas fraction between 0.0 and 1.0 (default 0.9)",
        type="float"
    )
    
    result.add_option(
        "-r", "--rscale", 
        default = 1.0,
        dest="rscale",
        help="length scale of the problem in parsec (default 1) ",
        type="float"
    )
    
    result.add_option(
        "--star_smoothing_fraction", 
        default = 0.001,
        dest="star_smoothing_fraction",
        help="smoothing length of the stars as a fraction of the length scale",
        type="float"
    )
    
    result.add_option(
        "--gas_smoothing_fraction", 
        default = 0.05,
        dest="gas_smoothing_fraction",
        help="smoothing length of the gas particles as a fraction of the length scale",
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
        "--interaction-timestep", 
        default = 0.01,
        dest="interaction_timestep",
        help="time between bridge interactions (0.01 nbody time)",
        type="float"
    )
    result.add_option(
        "-t", "--end-time", 
        default = 1,
        dest="endtime",
        help="end time of the simulation (in nbody time, default 1)",
        type="float"
    )
    result.add_option(
        "--ntimesteps", 
        default = 10,
        dest="ntimesteps",
        help="number of times to do reporting",
        type="int"
    )
    result.add_option(
        "--noplot", 
        dest="must_do_plot",
        default = True,
        help="do not show a plot and end as soon as possible",
        action="store_false"
    )
    result.add_option(
        "--allinoone", 
        dest="must_do_bridge",
        default = True,
        help="simulate the stars and gas with one sph code",
        action="store_false"
    )
    return result
    
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    if options.must_do_bridge:
        options = options.__dict__
        BridgeStarAndGasPlummerCode(**options)
    else:
        options = options.__dict__
        AllInOneStarAndGasPlummerCode(**options)
    
