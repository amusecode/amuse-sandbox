'''
This script simulates the effect of Bondi-Hoyle accretion when a star is emerged in a flow of ISM. It is possible to add a protoplanetary disk around the star
to see how this effects the accretion (onto the star). The computational domain comprises a cylinder. One can choose either the star to be fixed and the ISM 
to move through the cylinder, adding a new slice of ISM each timestep, or the star to be moving and the cylinder filled with ISM. 

If you have questions, do not hesitate to contact me: thomas.wijnen@astro.ru.nl
'''
import matplotlib 
matplotlib.use('Agg') #use Agg so python does not use the GUI for plotting. Must import pyplot after this command
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from amuse.community.gadget2_sphs.interface import Gadget2 as Gadget2visc
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.units import units, constants, nbody_system, generic_unit_system
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.units.quantities import VectorQuantity
from amuse.ext.sink import new_sink_particles
from amuse.datamodel import Particles, Grid
from amuse.io import write_set_to_file, read_set_from_file
from amuse.plot import pynbody_column_density_plot, scatter
from amuse.plot import plot as amuseplot
from amuse.units.optparse import OptionParser
from amuse.units.generic_unit_system import *
import numpy
import time as timing
import os
import glob


class CalculateExactSolutionIn1D(object):
    
    def __init__(self,     
                 cube_size=2. | generic_unit_system.length, #Size of the periodic box
                 rho1 = 1.0 | (generic_unit_system.mass/generic_unit_system.length**3), #speed of sound in the ISM
                 rho2 = 0.25 | (generic_unit_system.mass/generic_unit_system.length**3), #Velocity of the ISM with respect to the star
                 p1 = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))),
                 p2 = 0.1795 | (generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))),
                 gamma = 5./3):
        
        
        self.number_of_points = 1000
        
        self.rho1 = rho1#1.0 | density
        self.p1 = p1#1.0 | mass / (length * (time**2))
        self.u1 = 0.0 | speed
    
        self.rho5 = rho2# 0.25 | density
        self.p5 = p2#0.1795 | mass / (length * (time**2))
        self.u5 = 0.0 | speed
    
    
        self.gamma = gamma
    
        self.cube_length = cube_size
    
    def get_post_shock_pressure_p4(self, maximum_number_of_iterations = 20, maxium_allowable_relative_error = 1e-5):
        """solve for post-shock pressure by secant method"""
        p40 = self.p1
        p41 = self.p5
        p4  = p41
        
        f0  = self.calculate_p4_from_previous_value(p40)
        for x in range(maximum_number_of_iterations):
            f1 = self.calculate_p4_from_previous_value(p41)
            if (f1 == f0):
                return p4

            p4 = p41 - (p41 - p40) * f1 / (f1 - f0)

            error = abs (p4 - p41) / p41
            if (error < maxium_allowable_relative_error):
                return p4

            p40 = p41
            p41 = p4
            f0  = f1
        
        raise Exception("solution did not converge in less than {0!r} steps.".format(maximum_number_of_iterations))
        
    def get_post_shock_density_and_velocity_and_shock_speed(self, p4):
        z  = (p4 / self.p5 - 1.0)
        c5 = numpy.sqrt(self.gamma * self.p5 / self.rho5)
        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        gmfac1 = 0.5 * gm1 / self.gamma
        gmfac2 = 0.5 * gp1 / self.gamma

        fact = numpy.sqrt (1. + gmfac2 * z)

        u4 = c5 * z / (self.gamma * fact)
        rho4 = self.rho5 * (1.0 + gmfac2 * z) / (1.0 + gmfac1 * z)
        return u4, rho4, c5 * fact
        
        
    def calculate_p4_from_previous_value(self, p4):
        c1 = numpy.sqrt(self.gamma * self.p1 / self.rho1)
        c5 = numpy.sqrt(self.gamma * self.p5 / self.rho5)

        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        g2  = 2.0 * self.gamma
        
        z= (p4 / self.p5 - 1.0)
        fact = gm1 / g2 * (c5 / c1) * z / numpy.sqrt (1. + gp1 / g2 * z)
        fact = (1. - fact) ** (g2 / gm1)

        return self.p1 * fact - p4
        
    
    def get_solution_at_time(self, t):
        p4 = self.get_post_shock_pressure_p4()
        u4, rho4, w = self.get_post_shock_density_and_velocity_and_shock_speed(p4)
        
        #compute values at foot of rarefaction
        p3 = p4
        u3 = u4
        rho3 = self.rho1 * (p3 / self.p1)**(1. /self.gamma)
        
        c1 = numpy.sqrt (self.gamma * self.p1 / self.rho1)
        c3 = numpy.sqrt (self.gamma * p3 / rho3)
        
        xi = 0.0 | length#0.5 | length
        xr = -0.5*self.cube_length#1.0 | length
        xl = 0.5*self.cube_length#0.0 | length
        
        xsh = xi + w * t
        xcd = xi + u3 * t
        xft = xi + (u3 - c3) * t
        xhd = xi - c1 * t
        
        gm1 = self.gamma - 1.0
        gp1 = self.gamma + 1.0
        
        dx = (xr - xl) / (self.number_of_points - 1)
        x = xl + dx * numpy.arange(self.number_of_points)
        
        rho = VectorQuantity.zeros(self.number_of_points, density)
        p = VectorQuantity.zeros(self.number_of_points, mass / (length * time**2))
        u = VectorQuantity.zeros(self.number_of_points, speed)
        
        for i in range(self.number_of_points):
            if x[i] < xhd:
                rho[i] = self.rho1
                p[i]   = self.p1
                u[i]   = self.u1
            elif x[i] < xft:
                u[i]   = 2. / (self.gamma + 1.0) * (c1 + (x[i] - xi) / t)
                fact   = 1. - 0.5 * gm1 * u[i] / c1
                rho[i] = self.rho1 * fact ** (2. / gm1)
                p[i]   = self.p1 * fact ** (2. * self.gamma / gm1)
            elif x[i] < xcd:
                rho[i] = rho3
                p[i]   = p3
                u[i]   = u3
            elif x[i] < xsh:
                rho[i] = rho4
                p[i]   = p4
                u[i]   = u4
            else:
                rho[i] = self.rho5
                p[i]   = self.p5
                u[i]   = self.u5
                
        return x, rho,p,u



class SodTest(object):

    def __init__(self,     
                 total_N=5e4, #total number of ISM sph-particles when the cylinder is full
                 tend=0.1 | generic_unit_system.time, #End time of the simulation
                 cube_size=2. | generic_unit_system.length, #Size of the periodic box
                 dt = 0.001 | generic_unit_system.time, #timestep
                 rho1 = 1.0 | (generic_unit_system.mass/generic_unit_system.length**3), 
                 rho2 = 0.25 | (generic_unit_system.mass/generic_unit_system.length**3), 
                 p1 = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))),
                 p2 = 0.1795 | (generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))),
                 gamma = 1.4, 
                 verbose = False, #Control outputs to screen,
                 n_core = 4, #Number of cores used by the community codes
                 var_visc = True,
                 dirname = '/home/thomas/amuse-svn/examples/plots/sodtest',
                 code = "gadget"):
        
        self.total_N = total_N
        self.tend = tend
        self.dt = dt
        self.rho1 = rho1
        self.rho2 =rho2
        self.p1 = p1
        self.p2 = p2
        self.rho = min(self.rho1,self.rho2)
        self.ratio = max(self.rho1,self.rho2)/self.rho
        self.gamma = gamma
        self.u1 = (self.p1/self.rho1)/(gamma-1)
        self.u2 = (self.p2/self.rho2)/(gamma-1)
        self.code = code

        self.verbose = verbose
        self.var_visc = var_visc
        self.n_core = n_core #"Now can I get an encore, do you want more" - Jay-Z
        self.dirname = dirname
        
        self.converter = ConvertBetweenGenericAndSiUnits(1.| units.kpc, 1.e9 | units.MSun, constants.G)
        
        self.cube_length = cube_size
        self.cube_vol = self.cube_length**(3.)
        
        print "-------------------------------------"
        print "rho_1 = ", self.rho1.value_in((generic_unit_system.mass/generic_unit_system.length**3))
        print "rho_2 = ", self.rho2.value_in((generic_unit_system.mass/generic_unit_system.length**3))
        print "rho_1 in SI = ", self.converter.to_si(self.rho1)
        print "rho_2 in SI = ", self.converter.to_si(self.rho2)
        print "p_1 = ", self.p1.value_in((generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))))
        print "p_2 = ", self.p2.value_in((generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))))
        print "Internal energy 1 = ", self.u1.value_in((generic_unit_system.length**2/generic_unit_system.time**2))
        print "Internal energy 2 = ", self.u2.value_in((generic_unit_system.length**2/generic_unit_system.time**2))
        print 'Files will be saved in ', self.dirname 

    def setup(self):
        
        self.ism_code = self.create_codes()
        
        self.initialize_data()
     
    
    def create_codes(self):
        
        print "-------------------------------------"
        
        if self.code.find("gadget")>=0:
            if self.var_visc:
                ism_code = Gadget2visc(self.converter, number_of_workers=self.n_core,mode="periodic",redirection = 'none')#, debugger='gdb')
                print "Using Gadget2_viscosity"
            else:
                ism_code = Gadget2(self.converter, number_of_workers=self.n_core,mode="periodic")
                print "Using original Gadget2"
                
            ism_code.parameters.time_max = 10*self.tend





        
        elif self.code.find("fi")>=0:
            print "Using Fi"
            ism_code = Fi(self.converter, mode="periodic")#, redirection='none')#, debugger='gdb')
            ism_code.parameters.isothermal_flag=False
            ism_code.parameters.gamma=self.gamma
            #ism_code.parameters.integrate_entropy_flag=False
            ism_code.parameters.self_gravity_flag=False
            ism_code.parameters.periodic_box_size = self.cube_length
            ism_code.parameters.timestep = self.dt
            #print "Timestep is ", self.dt
            #ism_code.parameters.verbosity = 99
            #print ism_code.parameters
            
            

        else:
            print "No hydro code specified."
            return

        
        
        #ism_code.parameters.n_smooth = 200 #SPH codes adapt the smoothing length such that each particle will have 64 neighbours. Default is 50 in Amuse, 64 in Gadget.
        #ism_code.parameters.n_smooth_tol = 2./64. #fractional tolerance in number of SPH neighbours. Default is 2 neighbours in Gadget, fraction of 0.1 in Amuse   
        ism_code.parameters.artificial_viscosity_alpha = 0.1 #The artificial Bulk viscosity
        #ism_code.parameters.n_smooth_tol = 0.99999


        print "Number of cores used by community code(s): ", self.n_core
        print "Isothermal EoS = ", ism_code.parameters.isothermal_flag
        print "Artificial bulk viscosity is: ", ism_code.parameters.artificial_viscosity_alpha       
        print "Domain size = ", ism_code.parameters.periodic_box_size.value_in(units.AU), " AU"
        print "Fractional tolerance  = ", ism_code.parameters.n_smooth_tol
        print "Number of neighbours = ", ism_code.parameters.n_smooth
        
        return ism_code
        
    def initialize_data(self):
        
        self.n_particles = int((self.total_N / (self.ratio + 1)))

        self.low_particle_dens = self.n_particles/(0.5*self.cube_vol)
        
        self.particle_mass = self.rho/self.low_particle_dens
        
        print "-------------------------------------"       
        print " Cube length = ", self.cube_length.value_in(generic_unit_system.length)
        print " Cube volume = ", self.cube_vol.value_in((generic_unit_system.length**3.))
        print " Sph particle density in right half = ", self.low_particle_dens.value_in((generic_unit_system.length**(-3)))
        print " density of right half =  ", self.rho.value_in((generic_unit_system.mass/generic_unit_system.length**3))
        print " Sph particle mass = ", self.particle_mass.value_in(generic_unit_system.mass)
            
    def make_half(self, side):
    
        if side == 'left':
            n_particles = int(self.n_particles*self.ratio)
            
        elif side == 'right':
            n_particles = self.n_particles
        
    
        half = Particles(n_particles)
    

        n1D = numpy.ceil((n_particles)**(1./3) )
        x, y, z = numpy.mgrid[0. : 0.4999 : n1D*1j,
                                  0. : 0.999 : n1D*1j,
                                  0. : 0.999 : n1D*1j]
        half.x = x.flatten() * self.cube_length
        half.y = y.flatten() * self.cube_length
        half.z = z.flatten() * self.cube_length
    
    
        if side == 'left':
            half.u = self.u1
        elif side == 'right':
            half.x += 0.5*self.cube_length
            half.u = self.u2    
    
    
        #Particles only get no initial velocity
        half.vx = 0.0 | (generic_unit_system.length/generic_unit_system.time)
        half.vy = 0.0 | (generic_unit_system.length/generic_unit_system.time)
        half.vz = 0.0 | (generic_unit_system.length/generic_unit_system.time)   
    
        half.mass = self.particle_mass
        
        print "-------------------------------------"
        print "      "+str(side)+" hand side "
        print " x-range = ", half.x.min(), " - ", half.x.max()
        print " Number of particles = ", n_particles
        print " Particle mass = ", self.particle_mass.value_in(generic_unit_system.mass)
        print "Internal energy = ", half.u[0].value_in((generic_unit_system.length**2/generic_unit_system.time**2))
        print "Min internal energy  = ", self.converter.to_si(half.u.min())
        print "Max internal energy = ", self.converter.to_si(half.u.max())
        print "density  = ", ((self.particle_mass*n_particles)/(0.5*self.cube_vol)).value_in((generic_unit_system.mass/generic_unit_system.length**3))
    
        return half
    
    

    
    def evolve_model(self):
        
        left_particles = self.make_half('left')
        right_particles = self.make_half('right')
        right_particles.position -= [0.5, 0.5, 0.5] * self.cube_length
        left_particles.position -= [0.5, 0.5, 0.5] * self.cube_length


        
        self.leftside = self.ism_code.gas_particles.add_particles(left_particles)
        self.rightside = self.ism_code.gas_particles.add_particles(right_particles)
    
        time = 0.|generic_unit_system.time
        self.leading_zeros = len(str(int(self.tend/self.dt)))
        
        if self.var_visc:
            self.temp_instance = self.create_codes()
            
        self.exact = CalculateExactSolutionIn1D(cube_size=self.cube_length, #Size of the periodic box
                                                rho1 = self.rho1, #speed of sound in the ISM
                                                rho2 = self.rho2, #Velocity of the ISM with respect to the star
                                                p1 = self.p1,
                                                p2 = self.p2,
                                                gamma = self.gamma)
        
    
        while time <= self.tend:
                        
            if self.verbose: start = timing.time()
            print "======================================================="
            print "Evolving to time = ", time.value_in(generic_unit_system.time), " of ", self.tend.value_in(generic_unit_system.time)," seconds..."    
            self.ism_code.evolve_model(time)

            if self.verbose: print "This took ", (timing.time() - start), " s"

            if not self.var_visc: print "Artificial bulk viscosity is: ", self.ism_code.parameters.artificial_viscosity_alpha
            
            if not int(time/self.dt) % 2:
                self.make_plots(time)   

            time += self.dt 

            
        print "======================================================="
        self.ism_code.stop()
        if self.var_visc:
            self.temp_instance.stop()
        
        
    def plot_position(self, time):
        fig = pyplot.figure()
        pyplot.title(str(time.value_in(units.s))+" s")
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.leftside.x.value_in(units.m), self.leftside.y.value_in(units.m), self.leftside.z.value_in(units.m), c='r', marker='.')
        ax.scatter(self.rightside.x.value_in(units.m), self.rightside.y.value_in(units.m), self.rightside.z.value_in(units.m), c='b', marker='.')
        ax.set_xlabel('x axis (m)')
        ax.set_ylabel('y axis (m)')
        ax.set_zlabel('z axis (m)')

        ax.set_xlim3d(0, self.cube_length.value_in(units.m)) 
        ax.set_ylim3d(0, self.cube_length.value_in(units.m))
        ax.set_zlim3d(0, self.cube_length.value_in(units.m))            
    
        pyplot.savefig(self.dirname+'/position/sodtube'+str(int(time/self.dt)).zfill(self.leading_zeros)+'.eps')
        pyplot.close()

    def make_plots(self, time):  
        
        xpositions_exact, rho_exact, p_exact, u_exact = self.exact.get_solution_at_time(time)
        #xmin = -self.cube_length.value_in(generic_unit_system.length)
        xmax = self.cube_length.value_in(generic_unit_system.length)

        
        x = numpy.indices((int(xmax/0.0001 +1),))
        
        x = x.flatten()*0.0001-0.5*xmax        
        
        vx= 0.*x | (generic_unit_system.length/generic_unit_system.time) 
        vy=0.*vx
        vz=0.*vx
        
        
        
        x = generic_unit_system.length(x)
        y = 0.*x
        z = 0.*x
        
        if self.verbose: 
            print "    Getting hydro state at grid points..."
            start = timing.time()
        
        
        rho, rhovx, rhovy, rhovz, rhoe = [self.converter.to_generic(quantity) for quantity in self.ism_code.get_hydro_state_at_point(x, 
                                                                                y, 
                                                                                z, 
                                                                                vx, 
                                                                                vy, 
                                                                                vz)]
        
        
        if self.verbose: print "            This took ", (timing.time() - start), " s"
        
        
        pyplot.figure()
        pyplot.title(str(time.value_in(generic_unit_system.time))+" s")
        pyplot.xlabel('x')
        pyplot.ylabel('$\rho$')
        scatter(x, rho)
        amuseplot(xpositions_exact,rho_exact)
        pyplot.xlim(-0.5*xmax,0.5*xmax)         
        pyplot.ylim(0.,1.1*self.rho1.value_in((generic_unit_system.mass/generic_unit_system.length**3)))        
        pyplot.savefig(self.dirname+'/rho/rho'+str(int(time/self.dt)).zfill(self.leading_zeros)+'.eps')
        pyplot.close()
        
        
        pyplot.figure()
        pyplot.title(str(time.value_in(generic_unit_system.time))+" s")
        pyplot.xlabel('x')
        pyplot.ylabel('$v_{x}$')
        scatter(x, rhovx/rho)
        pyplot.xlim(-0.5*xmax,0.5*xmax)  
        pyplot.savefig(self.dirname+'/vx/vx'+str(int(time/self.dt)).zfill(self.leading_zeros)+'.eps')
        pyplot.close()
        
        
        #pressure_particles = self.ism_code.gas_particles.copy()
        visc_particles = self.ism_code.gas_particles.copy()
        '''    
        pressure_particles.mass *= self.converter.to_generic(pressure_particles.pressure).value_in((generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))))      
  
        self.temp_instance.gas_particles.add_particles(pressure_particles)
             
        print "    Calculating pressure at grid points..."
        start = timing.time()
        
        
        pres, presvx, presvy, presvz, prese = [self.converter.to_generic(quantity) for quantity in self.temp_instance.get_hydro_state_at_point(x, 
                                                                                     y, 
                                                                                     z, 
                                                                                     vx, 
                                                                                     vy, 
                                                                                     vz)]

        
        print "            This took ", (timing.time() - start), " s"
            
        self.temp_instance.gas_particles.remove_particles(pressure_particles)
        '''
        
        pres = (self.gamma-1)*(rhoe-0.5*((rhovx**2.)+(rhovy**2.)+(rhovz**2.))/rho)
        
        pyplot.figure()
        pyplot.title(str(time.value_in(generic_unit_system.time))+" s")
        scatter(x,pres)
        amuseplot(xpositions_exact,p_exact)
        pyplot.xlim(-0.5*xmax,0.5*xmax)
        pyplot.ylim(0.,1.1*self.p1.value_in((generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2)))))         
        pyplot.xlabel('x')
        pyplot.ylabel('P')
        pyplot.savefig(self.dirname+'/pres/pres'+str(int(time/self.dt)).zfill(self.leading_zeros)+'.eps')
        pyplot.close()
        
        
        if self.var_visc:
            visc_particles.mass *= visc_particles.alpha      
            
            self.temp_instance.gas_particles.add_particles(visc_particles)
             
            print "    Calculating viscosity at grid points..."
            start = timing.time()
            visc, viscvx, viscvy, viscvz, visce = [self.converter.to_generic(quantity) for quantity in self.temp_instance.get_hydro_state_at_point(x, 
                                                                                     y, 
                                                                                     z, 
                                                                                     vx, 
                                                                                     vy, 
                                                                                     vz)]
        
        
            print "            This took ", (timing.time() - start), " s"
            
            self.temp_instance.gas_particles.remove_particles(visc_particles)
        
            pyplot.figure()
            pyplot.title(str(time.value_in(generic_unit_system.time))+" s")
            scatter(x, visc/rho)#, linestyle ='dashed', color='black')
            pyplot.xlim(-0.5*xmax,0.5*xmax)
            pyplot.xlabel('x')
            pyplot.ylabel('Viscosity')
            pyplot.savefig(self.dirname+'/visc/visc'+str(int(time/self.dt)).zfill(self.leading_zeros)+'.eps')
            pyplot.close()
        
        
        
        
        

def new_option_parser():
    result = OptionParser()
    result.add_option("-n", 
                      dest="total_N", 
                      type="int", 
                      default = 5e5,
                      help="Total number of sph particles when cylinder is full")
    result.add_option("-t", 
                      dest="t_end", 
                      type="float", 
                      default = 0.1,
                      unit = generic_unit_system.time, 
                      help="end time [%unit]")
    result.add_option("-s", 
                      dest="size", 
                      type="float", 
                      default = 2.,
                      unit = generic_unit_system.length, 
                      help="Size of the periodic box [%unit]")
    result.add_option("-d", 
                      dest="dt_diag", 
                      type="float", 
                      default = 1./128.,
                      unit = generic_unit_system.time, 
                      help="diagnostic timestep [%unit]")
    result.add_option("--rho1", 
                      dest="rho1", 
                      type="float", 
                      default = 1.0,
                      unit = (generic_unit_system.mass/generic_unit_system.length**3),
                      help="Density of left hand")
    result.add_option("--rho2", 
                      dest="rho2", 
                      type="float", 
                      default = 0.25,
                      unit = (generic_unit_system.mass/generic_unit_system.length**3),
                      help="Density of right hand")
    result.add_option("--p1", 
                      dest="p1", 
                      type="float", 
                      default = 1.0,
                      unit = (generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))),
                      help="Pressure of left hand")
    result.add_option("--p2", 
                      dest="p2", 
                      type="float", 
                      default = 0.1795,
                      unit = (generic_unit_system.mass / (generic_unit_system.length * (generic_unit_system.time**2))),
                      help="Pressure of right hand") 
    result.add_option("-g", 
                      dest="gamma", 
                      type="float", 
                      default = 5./3,
                      help="Polytropic index (gamma)")
    result.add_option("--nover", 
                      dest="verbosity",
                      action='store_false', 
                      default = True,
                      help="Verbosity flag for output to the screen")
    result.add_option("--dir", 
                      dest="directory", 
                      type="string", 
                      default = '/home/thomas/amuse-svn/examples/plots/sodtest',
                      help="The directory where to save the plots and data")
    result.add_option("-c", 
                      dest="n_cores", 
                      type="int", 
                      default = 4,
                      help="Number of cores used by community codes")
    result.add_option("--novarvisc", 
                      dest="varvisc",  
                      action='store_false', 
                      default = True,
                      help="Use gadget with or without the variable viscosity")
    result.add_option("--code",
                      choices=["gadget", "fi"],
                      default="gadget",
                      dest="code",
                      help="Which code to use")
    return result


if __name__ in ("__main__","__plot__"):
    
    
    
    o, arguments  = new_option_parser().parse_args()   
    
    dirs = ['/pres', '/position', '/vx', '/rho', '/visc']
    
    for item in dirs: 
        
        if not os.path.exists(o.directory+item):
            print "Creating "+o.directory+item
            os.makedirs(o.directory+item)
        else:
        
            files = glob.glob(o.directory+item+"/*")
            if (len(files) > 0): print "Removing files from "+o.directory+item
            for f in files:
                os.remove(f)
        
    code = SodTest( total_N=o.total_N,
                       tend=o.t_end,
                       cube_size=o.size,
                       dt = o.dt_diag,
                       rho1 = o.rho1,
                       rho2 = o.rho2,
                       p1 = o.p1,
                       p2 = o.p2,
                       gamma = o.gamma,
                       verbose = o.verbosity,
                       dirname = o.directory,
                       n_core = o.n_cores,
                       var_visc = o.varvisc,
                       code = o.code)
    code.setup()
    code.evolve_model()
