from amuse.lab import *
from optparse import OptionParser
from amuse.io import store

def main(N=100, Mp=1.0, Ms=1.0, Rp=1.0, Rs=1.0, a=2.9, e=0.0, t_end=1, n_steps=10):
    t_end = t_end | units.day

    Mp = Mp |units.MSun
    Rp = Rp| units.RSun
    converterA=nbody_system.nbody_to_si(Mp, Rp)
    starA = new_plummer_gas_model(N, convert_nbody=converterA)
    starA.position += (0, 0, 0) | Rp.unit
    Ms = Ms |units.MSun
    Rs = Rs| units.RSun
    Ns = int(Ms/Mp * N)
    converterB=nbody_system.nbody_to_si(Ms, Rs)
    starB = new_plummer_gas_model(Ns, convert_nbody=converterB)
    starB.position += (a, 0, 0) | Rs.unit
    starA.h_smooth = 0|units.parsec
    starB.h_smooth = 0|units.parsec
    a = a | units.RSun
    r = a
    v_orb = (constants.G*(Mp+Ms)*(2./r - 1./a)).sqrt().value_in(units.kms)
    starB.velocity += (0, v_orb, 0) | units.kms
    stars = ParticlesSuperset([starA, starB])
    stars.move_to_center() 

    hydro = Gadget2(converterA)
    hydro.parameters.self_gravity_flag=1
    hydro.gas_particles.add_particles(stars)
    
    Etot_init = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    channel_from_hydro_to_framework = hydro.gas_particles.new_channel_to(starA)
    filename = "hydro.hdf5"
    time = 0.0 | t_end.unit
    write_set_to_file(hydro.particles, filename, 'hdf5')
    dt = t_end/float(n_steps)
    while time < t_end:
        time += dt

        hydro.evolve_model(time)
        channel_from_hydro_to_framework.copy()
        write_set_to_file(hydro.particles, filename, 'hdf5')

        Ekin = hydro.kinetic_energy 
        Epot = hydro.potential_energy
        Eth = hydro.thermal_energy
        Etot = Ekin + Epot + Eth
        print "T=", hydro.get_time(), "M=", hydro.gas_particles.mass.sum(), 
        print "E= ", Etot, "Q= ", (Ekin+Eth)/Epot, "dE=", (Etot_init-Etot)/Etot

    hydro.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of SPH particles [10]")
    result.add_option("-n", dest="n_steps", type="int",default = 10,
                      help="number of steps [10]")
    result.add_option("-t", dest="t_end", type="float", default = 1,
                      help="end time of the simulation [1] day")
    result.add_option("-M", dest="Mp", type="float", default = 1,
                      help="Mass of primary star [1] MSun")
    result.add_option("-m", dest="Ms", type="float", default = 1,
                      help="Mass of secondary star [1] MSun")
    result.add_option("-R", dest="Rp", type="float", default = 1,
                      help="Radiu of primary star [1] RSun")
    result.add_option("-r", dest="Rs", type="float", default = 1,
                      help="Radius of secondary star [1] RSun")
    result.add_option("-a", dest="a", type="float", default = 2.9,
                      help="orbital separation [1] Rsun")
    result.add_option("-e", dest="e", type="float", default = 0,
                      help="eccentricity [0]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

