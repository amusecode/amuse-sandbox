import numpy
import os
import os.path
import random
import sys
import pickle
from optparse import OptionParser
from time import clock

from amuse.community.ph4.interface import ph4 as GravityModule
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.couple import multiples

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import quantities

from amuse import datamodel
from amuse.datamodel import particle_attributes as pa
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody

import multiples_restart_functions as MRest

from amuse import io

def new_smalln():
    SMALLN.reset()
    return SMALLN

def get_coms_in_multiples(multiples_code):
    """ Returns the set of Centre of Mass particles in a Multiples module. """
    result = datamodel.Particles()
    for root, tree in multiples_code.root_to_tree.iteritems():
        root_particle = root.as_particle_in_set(multiples_code._inmemory_particles)
        result.add_particle(root)
    return result

def get_singles_in_multiples(multiples_code):
    """ Returns the set of individual particles in a Multiples module. """
    result = datamodel.Particles()
    for root, tree in multiples_code.root_to_tree.iteritems():
        root_particle = root.as_particle_in_set(multiples_code._inmemory_particles)
        #result.remove_particle(root)
        leaves = tree.get_leafs_subset()

        original_star = tree.particle

        dx = root_particle.x - original_star.x
        dy = root_particle.y - original_star.y
        dz = root_particle.z - original_star.z
        dvx = root_particle.vx - original_star.vx
        dvy = root_particle.vy - original_star.vy
        dvz = root_particle.vz - original_star.vz

        leaves_in_result = result.add_particles(leaves)
        leaves_in_result.x += dx
        leaves_in_result.y += dy
        leaves_in_result.z += dz
        leaves_in_result.vx += dvx
        leaves_in_result.vy += dvy
        leaves_in_result.vz += dvz
    return result

def print_log(pre, time, gravity, E0 = 0.0 | nbody_system.energy, cpu0 = 0.0):
    cpu = clock()
    N = len(gravity.particles)
    M = gravity.total_mass
    U = gravity.potential_energy
    T = gravity.kinetic_energy
    Etop = T + U
    Nmul, Nbin, Emul = gravity.get_total_multiple_energy()
    tmp1,tmp2,Emul2 = gravity.get_total_multiple_energy2()
    Etot = Etop + Emul
    Eext = gravity.multiples_external_tidal_correction
    Eint = gravity.multiples_internal_tidal_correction
    Eerr = gravity.multiples_integration_energy_error
    Edel = gravity.multiples_external_tidal_correction \
		+ gravity.multiples_internal_tidal_correction \
		+ gravity.multiples_integration_energy_error
    Ecor = Etot - Edel
    if E0 == 0 | nbody_system.energy: E0 = Ecor
    Rvir = -0.5*M*M/U
    Q = -T/U
    com = pa.center_of_mass(gravity.particles)
    comv = pa.center_of_mass_velocity(gravity.particles)
    dcen,rcore,rhocore = pa.densitycentre_coreradius_coredens(gravity.particles)
    cmx,cmy,cmz = dcen
    lagr,mf = pa.LagrangianRadii(gravity.particles, cm=dcen)  # no units!

    print ''
    print pre+"time=", time.number
    print pre+"CPU=", cpu - cpu0
    print pre+"Ntot=", N
    print pre+"mass=", M.number
    print pre+"Etot=", Etot.number
    print pre+"Etop=", Etop.number
    print pre+"Eext=", Eext.number
    print pre+"Eint=", Eint.number
    print pre+"Eerr=", Eerr.number
    print pre+"Edel=", Edel.number
    print pre+"Ecor=", Ecor.number
    print pre+"dE/E=", Ecor/E0 - 1
    print pre+"Rvir=", Rvir.number
    print pre+"Qvir=", Q
    cmx,cmy,cmz = com
    print pre+"cmpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number)
    cmx,cmy,cmz = comv
    print pre+"cmvel[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number)
    cmx,cmy,cmz = dcen
    print pre+"dcpos[3]= %.8f %.8f %.8f" % (cmx.number, cmy.number, cmz.number)
    print pre+"Rcore=", rcore.number
    print pre+"Mcore=", (rhocore*rcore**3).number	# fake...
    print pre+"Mlagr[9]=",
    for m in mf: print "%.4f" % m,
    print ''
    print pre+"Rlagr[9]=",
    for r in lagr.number: print "%.8f" % r,
    print ''
    kT = T/N
    #Nmul,Nbin,Emul = gravity.pretty_print_multiples(pre, kT, dcen)
    print pre+"Nmul=", Nmul
    print pre+"Nbin=", Nbin
    print pre+"Emul= %.5f" % Emul.number
    print pre+"Emul2= %.5f" % Emul2.number
    print pre+"Emul/kT= %.5f" % (Emul/kT)
    print pre+"Emul/E= %.5f" % (Emul/Etot)
    print ''

    sys.stdout.flush()
    return Ecor,cpu

SMALLN = None
def new_smalln():
    SMALLN.reset()
    return SMALLN

def init_smalln():
    global SMALLN
    #SMALLN = SmallN(redirection="none", debugger="xterm")
    SMALLN = SmallN(redirection="none")
    SMALLN.parameters.timestep_parameter = 0.1
    SMALLN.parameters.cm_index = 2001

def init_kepler(star1, star2):
    try:
        star1.mass.value_in(units.kg) # see if SI units, throw exception if not
        unit_converter \
            = nbody_system.nbody_to_si(star1.mass + star2.mass,
                                       (star2.position-star1.position).length())
    except Exception as ex:
        unit_converter = None

    kep = Kepler(unit_converter, redirection = "none")
    kep.initialize_code()

    return kep

def run_ph4(options, time=None, stars=None, mc_root_to_tree=None, randomize=True):
    infile = options.infile
    outfile = options.outfile
    restart_file = options.restart_file
    number_of_stars = options.N
    number_of_binaries = options.Nbin
    end_time = options.t_end | nbody_system.time
    delta_t = options.delta_t | nbody_system.time
    n_workers = options.n_workers
    use_gpu = options.use_gpu
    gpu_worker = options.gpu_worker
    salpeter = options.salpeter
    accuracy_parameter = options.accuracy_parameter
    softening_length = options.softening_length | nbody_system.length
    manage_encounters = options.manage_encounters
    random_seed = options.random_seed

    if randomize:
        if random_seed <= 0:
            numpy.random.seed()
            random_seed = numpy.random.randint(1, pow(2,31)-1)
        numpy.random.seed(random_seed)
        print "random seed =", random_seed

    if infile is not None: print "input file =", infile
    if restart_file is not None: print "restart file =", restart_file
    if restart_file is not None and infile is not None: print "restart file overrides input file"
    print "end_time =", end_time.number
    print "delta_t =", delta_t.number
    print "n_workers =", n_workers
    print "use_gpu =", use_gpu
    print "manage_encounters =", manage_encounters
    print "n =", number_of_stars
    print "nbin=", number_of_binaries
    print "\ninitializing the gravity module"
    sys.stdout.flush()

    init_smalln()

    # Note that there are actually three GPU options:
    #
    #	1. use the GPU code and allow GPU use (default)
    #	2. use the GPU code but disable GPU use (-g)
    #	3. use the non-GPU code (-G)

    if gpu_worker == 1:
        try:
            #gravity = GravityModule(number_of_workers = n_workers,
            #               redirection = "xterm")
            gravity = GravityModule(number_of_workers = n_workers,
                           redirection = "none", mode = "gpu")
        except Exception as ex:
            gravity = GravityModule(number_of_workers = n_workers,
                           redirection = "none")
    else:
        gravity = GravityModule(number_of_workers = n_workers,
                       redirection = "none")

    gravity.initialize_code()
    gravity.parameters.set_defaults()

    if softening_length < 0.0 | nbody_system.length:

        # Use ~interparticle spacing.  Assuming standard units here.  TODO

        eps2 = 0.25*(float(number_of_stars))**(-0.666667) \
            | nbody_system.length**2
    else:
        eps2 = softening_length*softening_length
    print 'softening length =', eps2.sqrt()

    gravity.parameters.timestep_parameter = accuracy_parameter
    gravity.parameters.epsilon_squared = eps2
    gravity.parameters.use_gpu = use_gpu    

    kep = Kepler(redirection = "none")
    kep.initialize_code()

    multiples_code = None

    #-----------------------------------------------------------------

    if (restart_file is None or not os.path.exists(restart_file+".stars.hdf5")) and infile is None and stars is None:

        print "making a Plummer model"
        stars = new_plummer_model(number_of_stars)

        id = numpy.arange(number_of_stars)
        stars.id = id+1

        print "setting particle masses and radii"
        if salpeter == 0:
            print 'equal masses'
            total_mass = 1.0 | nbody_system.mass
            scaled_mass = total_mass / number_of_stars
        else:
            print 'salpeter mass function'
            scaled_mass = new_salpeter_mass_distribution_nbody(number_of_stars)
        stars.mass = scaled_mass

        print "centering stars"
        stars.move_to_center()
        print "scaling stars to virial equilibrium"
        stars.scale_to_standard(smoothing_length_squared
                                    = gravity.parameters.epsilon_squared)
        time = 0.0 | nbody_system.time

        total_mass = stars.mass.sum()
        ke = pa.kinetic_energy(stars)
        kT = ke/(1.5*number_of_stars)

        if number_of_binaries > 0:

            # Turn selected stars into binary components.
            # Only tested for equal-mass case.

            added_mass = 0.0 | nbody_system.mass

            # Work with energies rather than semimajor axes.

            Emin = 10*kT
            Emax = 20*kT
            ecc = 0.1

            id_count = number_of_stars
            nbin = 0
            for i in range(0, number_of_stars,
                           number_of_stars/number_of_binaries):

                # Star i is CM, becomes component, add other star at end.

                nbin += 1

                mass = stars[i].mass
                #new_mass = numpy.random.uniform()*mass	# uniform q?
                new_mass = mass	# uniform q?
                mbin = mass + new_mass
                fac = new_mass/mbin
                E = Emin + numpy.random.uniform()*(Emax-Emin)
                a = 0.5*nbody_system.G*mass*new_mass/E

                kep.initialize_from_elements(mbin, a, ecc)
                dr = quantities.AdaptingVectorQuantity()
                dr.extend(kep.get_separation_vector())
                dv = quantities.AdaptingVectorQuantity()
                dv.extend(kep.get_velocity_vector())

                newstar = datamodel.Particles(1)
                newstar.mass = new_mass
                newstar.position = stars[i].position + (1-fac)*dr
                newstar.velocity = stars[i].velocity + (1-fac)*dv
                # stars[i].mass = mass
                stars[i].position = stars[i].position - fac*dr
                stars[i].velocity = stars[i].velocity - fac*dv

                id_count += 1
                newstar.id = id_count
                stars.add_particles(newstar)
                added_mass += new_mass

                if nbin >= number_of_binaries: break

            kep.stop()

            print 'created', nbin, 'binaries'
            sys.stdout.flush()

            stars.mass = stars.mass * total_mass/(total_mass+added_mass)
            number_of_stars += nbin
            Xtra=numpy.zeros(2)
        print "recentering stars"
        stars.move_to_center()
        sys.stdout.flush()

        stars.savepoint(time)
        
        print ''
        print "adding particles"
        # print stars
        sys.stdout.flush()
        gravity.particles.add_particles(stars)
        gravity.commit_particles()
    else:
        print "Restart detected.  Loading parameters from restart."
        new_end = options.t_end
        stars, time, multiples_code, Xtra = MRest.read_state_from_file(restart_file, gravity, new_smalln, kep)
        options.t_end = new_end
    
    total_mass = stars.mass.sum()
    ke = pa.kinetic_energy(stars)
    kT = ke/(1.5*number_of_stars)

    # Set dynamical radii (assuming virial equilibrium and standard
    # units).  Note that this choice should be refined, and updated
    # as the system evolves.  Probably the choice of radius should be
    # made entirely in the multiples module.  TODO.  In these units,
    # M = 1 and <v^2> = 0.5, so the mean 90-degree turnaround impact
    # parameter is
    #
    #		b_90 = G (m_1+m_2) / vrel^2
    #		     = 2 <m> / 2<v^2>
    #		     = 2 / N			for equal masses
    #
    # Taking r_i = m_i / 2<v^2> = m_i in virial equilibrium means
    # that, approximately, "contact" means a 90-degree deflection (r_1
    # + r_2 = b_90).  A more conservative choice with r_i less than
    # this value will isolates encounters better, but also place more
    # load on the large-N dynamical module.

    stars.radius = stars.mass.number | nbody_system.length

    # print "IDs:", stars.id.number

    print ''
    print "number_of_stars =", number_of_stars
    print "evolving to time =", end_time.number, \
          "in steps of", delta_t.number
    sys.stdout.flush()

    # Channel to copy values from the code to the set in memory.
    channel = gravity.particles.new_channel_to(stars)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()

    # -----------------------------------------------------------------
    # Create the coupled code and integrate the system to the desired
    # time, managing interactions internally.

    kep = init_kepler(stars[0], stars[1])
    if not multiples_code:
        multiples_code = multiples.Multiples(gravity, new_smalln, kep)

        multiples_code.neighbor_distance_factor = 1.0
        multiples_code.neighbor_veto = True
        #multiples_code.neighbor_distance_factor = 2.0
        #multiples_code.neighbor_veto = True
        multiples_code.retain_binary_apocenter = False

    print ''
    print 'multiples_code.initial_scale_factor =', \
        multiples_code.initial_scale_factor
    print 'multiples_code.neighbor_distance_factor =', \
        multiples_code.neighbor_distance_factor
    print 'multiples_code.neighbor_veto =', \
        multiples_code.neighbor_veto
    print 'multiples_code.final_scale_factor =', \
        multiples_code.final_scale_factor
    print 'multiples_code.initial_scatter_factor =', \
        multiples_code.initial_scatter_factor
    print 'multiples_code.final_scatter_factor =', \
        multiples_code.final_scatter_factor
    print 'multiples_code.retain_binary_apocenter =', \
        multiples_code.retain_binary_apocenter

#    if mc_root_to_tree is not None:
#        multiples_code.root_to_tree = mc_root_to_tree
#        print 'multiples code re-loaded with binary trees snapshot'

    pre = "%%% "
    E0,cpu0 = print_log(pre, time, multiples_code)
    
    while time < end_time:

        time += delta_t
        multiples_code.evolve_model(time)

        # Copy values from the module to the set in memory.

        channel.copy()

        # Copy the index (ID) as used in the module to the id field in
        # memory.  The index is not copied by default, as different
        # codes may have different indices for the same particle and
        # we don't want to overwrite silently.

        channel.copy_attribute("index_in_code", "id")

        print_log(pre, time, multiples_code, E0, cpu0)
        stars.savepoint(time)
        MRest.write_state_to_file(time, stars, gravity, multiples_code, options.restart_file, Xtra, backup=1)
        sys.stdout.flush()

    #-----------------------------------------------------------------

    if not outfile is None:

        # Write data to a file.

        f = open(outfile, 'w')

        #--------------------------------------------------
        # Need to save top-level stellar data and parameters.
        # Need to save multiple data and parameters.

        f.write('%.15g\n'% time.number)
        for s in multiples_code.stars: write_star(s, f)

        #--------------------------------------------------

        f.close()
        print 'wrote file', outfile

    print ''
    gravity.stop()


def write_star(s, f):
    x,y,z = s.position.number
    vx,vy,vz = s.velocity.number
    f.write('%d %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n' \
		%(s.id, s.mass.number, x, y, z, vx, vy, vz))

if __name__ == '__main__':

    print '\ncommand line:',
    for a in sys.argv: print a,
    print '\n'

    parser = OptionParser()
    parser.set_defaults(use_gpu=False)
    parser.add_option("-a", "--accuracy-parameter", dest="accuracy_parameter", 
                      default=0.1, type="float", 
                      help="Accuracy parameter for the top-level N-Body code.")
    parser.add_option("-b", "--num-binaries", dest="Nbin", default=0, type="int",
                      help="Number of binary stars to create at initialization.")
    parser.add_option("-t", "--tend", dest="t_end", default=5.0, type="float",
                      help="N-Body time at which to stop the simulation.")
    parser.add_option("-d", "--deltat", dest="delta_t", default=1.0,
                      type="float", help="Time interval for output.")
    parser.add_option("-f", "--infile", dest="infile", default=None,
                      type="string", help="Initial conditions file.")
    parser.add_option("-F", "--outfile", dest="outfile", default=None,
                      type="string", help="Final conditions file.")
    parser.add_option("-R", "--restartfile", dest="restart_file", default=None,
                      type="string", help="File for holding restart data.")
    parser.add_option("-n", "--num-singles", dest="N", default=100, type="int",
                      help="Number of single stars to create at initialization.")
    parser.add_option("-w", "--num-workers", dest="n_workers", default=2,
                      type="int", help="Number of worker processes to create.")
    parser.add_option("-g", "--use-gpu", dest="use_gpu", action="store_const",
                      const=1, help="Enable GPU for computation.")
    parser.add_option("-G", "--no-gpu", dest="use_gpu", action="store_const",
                      const=0, help="Disable GPU for computation.")
    parser.add_option("-W", "--num-gpu-workers", dest="gpu_worker", type="int",
                      default=1, help="Number of GPU workers.")
    parser.add_option("-S", "--salpeter", dest="salpeter", default=0,
                      action="store_const", const=1, 
                      help="Use a Salpeter mass function.")
    parser.add_option("-e", "--softening-length", dest="softening_length",
                      type="float", default=0, 
                      help="Softening length for top-level code.")
    parser.add_option("-s", "--random-seed", type="int", default=-1,
                      help="Random number seed.")
    parser.add_option("-c", "--manage-encounters", type="int", default=1,
                      help="Turn on/off encounter management.")

    (options, args) = parser.parse_args()

    assert is_mpd_running()

    run_ph4(options)

